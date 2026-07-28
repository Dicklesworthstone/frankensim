//! G0/G3/G4/G5 evidence for the solve orchestration driver
//! (bead frankensim-extreal-program-f85xj.6.5, slice 1).
//!
//! The battery drives the library seam directly: fixture project, real
//! import into an in-memory ledger, then the staged solve engine with a
//! scripted clock and caller-owned cancellation gate.

use fs_cli::{
    CompletedStage, GeometryImportLimits, GeometryImportRun, MAX_SOLVE_VISIBLE_OP_IDS,
    RawGeometryLibrary, SOLVE_DRIVER_VERSION, SolveCancellationPlan, SolveDriverState,
    SolveEvidencePhase, SolveRefusal, SolveRunId, SolveRunStatus, SolveStage,
    import_project_geometry, resume_solve, resume_solve_with_cancellation_plan, run_solve,
    run_solve_with_cancellation_plan,
};
use fs_exec::solver::LegacySnapshotV1Adapter;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::{NamedFaceGroup, quarantine::import_mesh};
use fs_ledger::{
    CONTROLLED_ARTIFACT_TILE_LEN, CONTROLLED_OP_FIELD_TILE_LEN, EdgeRole, ExtensionTable,
    FiveExplicits, Ledger, MAX_OP_FIELD_BYTES, OpOutcome, STORAGE_CHUNK_LEN, hash_bytes,
};
use fs_project::{
    Budgets, ConsequenceClass, Cooling, DecisionGate, DecodedProject, EntityDecl, Envelope,
    GeometryArtifact, GeometryAssignment, HalfSpaceSide, MeshSelector, Metadata, OutputRequest,
    PowerDissipation, ProjectSpec, RequirementDirection, RequirementSeverity, RequirementSource,
    RequirementSourceKind, SafetyFactorPolicy, Seeds, SolverSettings, ThermalLimit, UnitsDoctrine,
    Versions, print_sexpr,
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

fn tetra_stl() -> Vec<u8> {
    let p0 = [0.0, 0.0, 0.0];
    let p1 = [1.0, 0.0, 0.0];
    let p2 = [0.0, 1.0, 0.0];
    let p3 = [0.0, 0.0, 1.0];
    let mut stl = String::from("solid enclosure\n");
    stl.push_str(&facet(p0, p2, p1));
    stl.push_str(&facet(p0, p1, p3));
    stl.push_str(&facet(p0, p3, p2));
    stl.push_str(&facet(p1, p2, p3));
    stl.push_str("endsolid enclosure\n");
    stl.into_bytes()
}

fn storage_row_spanning_tetra_stl() -> Vec<u8> {
    let base = tetra_stl();
    let split = base
        .windows(b"facet normal".len())
        .position(|window| window == b"facet normal")
        .expect("fixture contains a facet");
    let mut padded = Vec::with_capacity(base.len() + STORAGE_CHUNK_LEN + 17);
    padded.extend_from_slice(&base[..split]);
    padded.resize(padded.len() + STORAGE_CHUNK_LEN + 17, b'\n');
    padded.extend_from_slice(&base[split..]);
    padded
}

fn project_for_receipt(seed_root: u64, source_hash: u64, parser_version: &str) -> ProjectSpec {
    let kelvin = |value| QtyAny::new(value, fs_project::spec::dims::TEMPERATURE);
    let watts = |value| QtyAny::new(value, fs_project::spec::dims::POWER);
    ProjectSpec {
        metadata: Some(Metadata {
            name: "solve-reference".to_string(),
            created: "2026-07-26".to_string(),
            context_of_use: "solve orchestration conformance".to_string(),
            intended_decision: "exercise the staged solve driver".to_string(),
            decision_gate: DecisionGate::ScopingEstimate,
            consequence: ConsequenceClass::Advisory,
        }),
        versions: Some(Versions {
            schema: fs_project::FSIM_VERSION,
            constellation: "00".repeat(32),
            workspace: "11".repeat(20),
        }),
        seeds: Some(Seeds { root: seed_root }),
        budgets: Some(Budgets {
            solve_time: QtyAny::new(60.0, fs_project::spec::dims::TIME),
            memory_bytes: 64 * 1024 * 1024,
            accuracy_rel: 0.01,
        }),
        capabilities: Some(vec!["thermal.conduction-solve".to_string()]),
        units: Some(UnitsDoctrine {
            storage: "si-base".to_string(),
            display: "engineering".to_string(),
        }),
        geometry: Some(vec![GeometryArtifact {
            role: "enclosure".to_string(),
            format: "stl".to_string(),
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
        perfect_contacts: None,
        power: Some(vec![PowerDissipation {
            region: "air".to_string(),
            watts: watts(5.0),
            duty: 1.0,
        }]),
        cooling: Some(Cooling {
            fans: Vec::new(),
            vents: Vec::new(),
            leakage: watts(0.0),
            airflow_leakage: None,
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
                document: "solve-fixture".to_string(),
                version: "1".to_string(),
                locator: "temperature-max".to_string(),
            },
            safety_factor: SafetyFactorPolicy {
                factor: 1.0,
                source: RequirementSource {
                    kind: RequirementSourceKind::UserDeclaration,
                    document: "solve-fixture-margin-policy".to_string(),
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

fn fixture_project(seed_root: u64, bytes: &[u8]) -> ProjectSpec {
    let receipt = import_mesh(bytes, "stl")
        .expect("fixture parses")
        .source_receipt;
    project_for_receipt(seed_root, receipt.source_hash, receipt.parser_version)
}

fn decode(spec: &ProjectSpec) -> DecodedProject {
    let source = print_sexpr(spec).expect("fixture renders canonically");
    let decoded = fs_project::parse_sexpr(&source).expect("fixture re-parses strictly");
    assert!(decoded.findings().is_empty(), "fixture validates cleanly");
    decoded
}

/// Import the fixture geometry into the ledger so the solve prefix has
/// retained evidence to verify.
fn import_fixture(ledger: &Ledger, spec: &ProjectSpec, bytes: Vec<u8>) -> GeometryImportRun {
    let artifact = &spec.geometry.as_ref().expect("geometry")[0];
    let mut raw = RawGeometryLibrary::new();
    assert!(!raw.insert_mesh(
        artifact,
        "fixtures/enclosure.stl",
        bytes,
        "m",
        0,
        vec![NamedFaceGroup {
            name: "fixture-all-faces".to_string(),
            // Deliberately non-sorted: solve must preserve writer order while
            // checking duplicates without a post-parse sort.
            faces: vec![3, 1, 2, 0],
        }],
    ));
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        import_project_geometry(spec, &raw, ledger, GeometryImportLimits::DEFAULT, cx)
            .expect("fixture imports")
    })
}

/// Reproduce a successful import candidate in an isolated ledger while
/// allowing one exact IR or summary-byte mutation.
fn reproduce_import_candidate(
    source: &Ledger,
    imported: &GeometryImportRun,
    target: &Ledger,
    ir: &str,
    summary: &str,
) -> i64 {
    reproduce_import_candidate_with_overrides(source, imported, target, ir, summary, None, None)
}

fn reproduce_import_candidate_with_mesh(
    source: &Ledger,
    imported: &GeometryImportRun,
    target: &Ledger,
    ir: &str,
    summary: &str,
    promoted_mesh_override: Option<&[u8]>,
) -> i64 {
    reproduce_import_candidate_with_overrides(
        source,
        imported,
        target,
        ir,
        summary,
        promoted_mesh_override,
        None,
    )
}

fn reproduce_import_candidate_with_overrides(
    source: &Ledger,
    imported: &GeometryImportRun,
    target: &Ledger,
    ir: &str,
    summary: &str,
    promoted_mesh_override: Option<&[u8]>,
    promotion_receipt_override: Option<&[u8]>,
) -> i64 {
    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let op = target
        .begin_op(
            None,
            ir,
            &FiveExplicits {
                seed: &source_row.seed,
                versions: &source_row.versions,
                budget: &source_row.budget,
                capability: &source_row.capability,
            },
            source_row.t_start,
        )
        .expect("reproduced import op");
    for entry in &imported.artifacts {
        for (hash, kind, role) in [
            (entry.raw_source, "geometry-source", EdgeRole::In),
            (
                entry.promotion_receipt,
                "geometry-import-receipt",
                EdgeRole::Out,
            ),
            (entry.promoted_mesh, "geometry-mesh-ply", EdgeRole::Out),
            (
                entry.assignment_report,
                "geometry-assignment-report",
                EdgeRole::Out,
            ),
        ] {
            let override_bytes = match kind {
                "geometry-mesh-ply" => promoted_mesh_override,
                "geometry-import-receipt" => promotion_receipt_override,
                _ => None,
            };
            let bytes = override_bytes.map_or_else(
                || {
                    source
                        .get_artifact(&hash)
                        .expect("read source artifact")
                        .expect("source artifact exists")
                },
                <[u8]>::to_vec,
            );
            let copied = target
                .put_artifact(kind, &bytes, None)
                .expect("copy typed import artifact");
            if override_bytes.is_some() {
                assert_eq!(copied.hash, hash_bytes(&bytes));
            } else {
                assert_eq!(copied.hash, hash);
            }
            target.link(op, &copied.hash, role).expect("typed edge");
        }
        let promotion_receipt = String::from_utf8(
            source
                .get_artifact(&entry.promotion_receipt)
                .expect("read source promotion receipt")
                .expect("source promotion receipt exists"),
        )
        .expect("source promotion receipt UTF-8");
        target
            .put_extension(
                ExtensionTable::Imports,
                &entry.import_record,
                &promotion_receipt,
            )
            .expect("copy exact import extension row");
    }
    let summary = target
        .put_artifact("geometry-import-run-receipt", summary.as_bytes(), None)
        .expect("reproduced summary");
    target
        .link(op, &summary.hash, EdgeRole::Out)
        .expect("summary edge");
    target
        .finish_op(
            op,
            OpOutcome::Ok,
            None,
            source_row.t_end.expect("source import finish clock"),
        )
        .expect("finish reproduced import");
    op
}

fn replace_json_object_with_empty(input: &str, marker: &str) -> String {
    let start = input.find(marker).expect("JSON marker") + marker.len();
    assert_eq!(input.as_bytes().get(start), Some(&b'{'));
    let mut depth = 0usize;
    let mut in_string = false;
    let mut escaped = false;
    let mut end = None;
    for (offset, byte) in input.as_bytes()[start..].iter().copied().enumerate() {
        if in_string {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                in_string = false;
            }
            continue;
        }
        match byte {
            b'"' => in_string = true,
            b'{' => depth += 1,
            b'}' => {
                depth = depth.checked_sub(1).expect("balanced JSON object");
                if depth == 0 {
                    end = Some(start + offset + 1);
                    break;
                }
            }
            _ => {}
        }
    }
    let end = end.expect("complete JSON object");
    format!("{}{{}}{}", &input[..start], &input[end..])
}

fn replace_json_unsigned(input: &str, marker: &str, replacement: usize) -> String {
    let start = input.find(marker).expect("JSON unsigned marker") + marker.len();
    let digits = input.as_bytes()[start..]
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    assert!(digits > 0, "marker must precede an unsigned JSON integer");
    format!(
        "{}{replacement}{}",
        &input[..start],
        &input[start + digits..]
    )
}

/// A benign clock: strictly increasing, one millisecond per call.
fn benign_clock() -> impl FnMut() -> f64 {
    let mut calls = 0u64;
    move || {
        calls += 1;
        #[allow(clippy::cast_precision_loss)]
        {
            calls as f64 * 0.001
        }
    }
}

fn run_to_gap(ledger: &Ledger, decoded: &DecodedProject) -> (SolveRefusal, Vec<String>) {
    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve(ledger, &gate, &mut clock, decoded, &mut progress)
        .expect_err("slice 1 refuses at the first stage gap");
    (refusal, progress)
}

fn run_to_one_stage_prefix(
    ledger: &Ledger,
    decoded: &DecodedProject,
) -> (SolveRefusal, Vec<String>) {
    let gate = CancelGate::new_clock_free();
    let gate_ref = &gate;
    let mut calls = 0u64;
    let mut clock = move || {
        calls += 1;
        if calls == 2 {
            gate_ref.request();
        }
        #[allow(clippy::cast_precision_loss)]
        {
            calls as f64 * 0.001
        }
    };
    let mut progress = Vec::new();
    let refusal = run_solve(ledger, &gate, &mut clock, decoded, &mut progress)
        .expect_err("caller gate stops after the first durable stage");
    assert_eq!(refusal.code, "cli-solve-cancelled");
    (refusal, progress)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SolvePublicationCounts {
    ops: u64,
    artifacts: u64,
    artifact_chunks: u64,
    edges: u64,
    artifact_output_seals: u64,
    op_artifact_edge_seals: u64,
    op_content_identities: u64,
    imports: u64,
}

fn solve_publication_counts(ledger: &Ledger) -> SolvePublicationCounts {
    SolvePublicationCounts {
        ops: ledger.table_count("ops").expect("op count"),
        artifacts: ledger.table_count("artifacts").expect("artifact count"),
        artifact_chunks: ledger
            .table_count("artifact_chunks")
            .expect("artifact chunk count"),
        edges: ledger.table_count("edges").expect("edge count"),
        artifact_output_seals: ledger
            .table_count("artifact_output_seals")
            .expect("artifact-output-seal count"),
        op_artifact_edge_seals: ledger
            .table_count("op_artifact_edge_seals")
            .expect("op-edge-seal count"),
        op_content_identities: ledger
            .table_count("op_content_identities")
            .expect("op-content-identity count"),
        imports: ledger.table_count("imports").expect("import count"),
    }
}

#[test]
fn g0_run_identity_is_deterministic_and_input_sensitive() {
    assert_eq!(
        SOLVE_DRIVER_VERSION, 2,
        "authority-semantic changes must deliberately advance this identity-bearing version"
    );

    let bytes = tetra_stl();
    let base = decode(&fixture_project(7, &bytes));
    let again = decode(&fixture_project(7, &bytes));
    assert_eq!(
        SolveRunId::derive(&base).to_hex(),
        SolveRunId::derive(&again).to_hex(),
        "identical projects derive identical run ids"
    );

    let seed_moved = decode(&fixture_project(8, &bytes));
    assert_ne!(
        SolveRunId::derive(&base).to_hex(),
        SolveRunId::derive(&seed_moved).to_hex(),
        "the RNG root seed is identity-bearing"
    );

    let mut workspace_moved_spec = fixture_project(7, &bytes);
    workspace_moved_spec
        .versions
        .as_mut()
        .expect("versions")
        .workspace = "22".repeat(20);
    let workspace_moved = decode(&workspace_moved_spec);
    assert_ne!(
        SolveRunId::derive(&base).to_hex(),
        SolveRunId::derive(&workspace_moved).to_hex(),
        "the declared workspace version is identity-bearing"
    );
}

#[test]
fn g0_solve_executes_the_real_prefix_then_refuses_at_the_first_gap() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let ops_after_import = ledger.table_count("ops").expect("count");

    let (refusal, progress) = run_to_gap(&ledger, &decoded);
    assert_eq!(refusal.code, "cli-solve-stage-gap");
    assert_eq!(refusal.stage, Some("material-resolve"));
    assert_eq!(refusal.dependency, Some("frankensim-hp7tb"));
    assert!(refusal.recorded_op.is_some(), "the gap refusal is ledgered");
    let run = refusal.run.clone().expect("run id derived");
    assert_eq!(run, SolveRunId::derive(&decoded).to_hex());

    // Two completed stage ops plus one recorded refusal op landed.
    let ops_after_solve = ledger.table_count("ops").expect("count");
    assert_eq!(ops_after_solve, ops_after_import + 3);

    // Both real stages reported progress.
    assert!(progress.iter().any(|line| line.contains("import-verify")));
    assert!(progress.iter().any(|line| line.contains("\"assign\"")));

    // Every solve op carries the run identity as its session.
    let run_id = SolveRunId::parse_hex(&run).expect("hex");
    let ids = ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .expect("ops");
    let solve_ops: Vec<i64> = ids
        .into_iter()
        .filter(|id| {
            ledger
                .op(*id)
                .expect("op row")
                .is_some_and(|row| row.session.as_deref() == Some(run_id.as_bytes().as_slice()))
        })
        .collect();
    assert_eq!(solve_ops.len(), 3, "two stages plus the recorded refusal");
}

#[test]
fn g3_solve_without_retained_import_refuses_with_evidence_diagnosis() {
    let bytes = tetra_stl();
    let decoded = decode(&fixture_project(7, &bytes));
    let ledger = Ledger::open(":memory:").expect("ledger");

    let (refusal, _) = run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
    assert_eq!(refusal.stage, Some("import-verify"));
    assert!(refusal.what.contains("no completed geometry import"));
    assert!(refusal.fix.contains("frankensim import"));
    assert!(
        refusal.recorded_op.is_some(),
        "evidence refusal is ledgered"
    );
}

fn run_to_gap_expect_code(
    ledger: &Ledger,
    decoded: &DecodedProject,
    code: &str,
) -> (SolveRefusal, Vec<String>) {
    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal =
        run_solve(ledger, &gate, &mut clock, decoded, &mut progress).expect_err("refusal expected");
    assert_eq!(refusal.code, code, "{refusal:?}");
    (refusal, progress)
}

#[test]
fn g3_an_import_for_a_different_project_does_not_satisfy_verification() {
    let bytes = tetra_stl();
    let imported_spec = fixture_project(7, &bytes);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &imported_spec, bytes.clone());

    // Same geometry, different seed: a different project hash, so the
    // retained import must not satisfy this run.
    let other = decode(&fixture_project(8, &bytes));
    let (refusal, _) = run_to_gap_expect_code(&ledger, &other, "cli-solve-import-evidence");
    assert!(refusal.what.contains("no completed geometry import"));
}

#[test]
fn g4_a_precancelled_solve_publishes_nothing() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let ops_before = ledger.table_count("ops").expect("count");
    let artifacts_before = ledger.table_count("artifacts").expect("count");

    let gate = CancelGate::new_clock_free();
    gate.request();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve(&ledger, &gate, &mut clock, &decoded, &mut progress)
        .expect_err("pre-cancelled run refuses");
    assert_eq!(refusal.code, "cli-solve-cancelled");
    assert!(refusal.recorded_op.is_none());

    assert_eq!(ledger.table_count("ops").expect("count"), ops_before);
    assert_eq!(
        ledger.table_count("artifacts").expect("count"),
        artifacts_before
    );
}

#[test]
fn g4_in_stage_evidence_cancellation_preserves_atomic_prefixes_at_every_owned_phase() {
    let phases = [
        (SolveEvidencePhase::ProjectValidation, None, 0),
        (SolveEvidencePhase::ProjectValidation, None, 1),
        (SolveEvidencePhase::ProjectIdentityDerive, None, 0),
        (SolveEvidencePhase::ProjectIdentityDerive, None, 1),
        (SolveEvidencePhase::FiveExplicitsRender, None, 0),
        (SolveEvidencePhase::FiveExplicitsRender, None, 1),
        (SolveEvidencePhase::CanonicalProjectRender, None, 0),
        (SolveEvidencePhase::CanonicalProjectRender, None, 1),
        (SolveEvidencePhase::ProjectIdentityDerive, Some(0), 0),
        (SolveEvidencePhase::ProjectIdentityDerive, Some(0), 1),
        (SolveEvidencePhase::VisibleOpPage, Some(0), 0),
        (SolveEvidencePhase::VisibleOpPage, Some(0), 1),
        (SolveEvidencePhase::VisibleOpPage, Some(0), u64::MAX),
        (SolveEvidencePhase::CandidateOpRowRead, Some(0), 0),
        (SolveEvidencePhase::CandidateOpRowRead, Some(0), u64::MAX),
        (SolveEvidencePhase::CandidateOpTextConversion, Some(0), 0),
        (SolveEvidencePhase::CandidateOpTextConversion, Some(0), 1),
        (SolveEvidencePhase::ImportIrParse, None, 1),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 0),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 1),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 2),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 0),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 1),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 17),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), 0),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), 1),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), u64::MAX),
        (SolveEvidencePhase::OperationContentIdentity, Some(0), 0),
        (
            SolveEvidencePhase::OperationContentIdentity,
            Some(0),
            u64::MAX,
        ),
        (SolveEvidencePhase::EdgePageRead, Some(0), 0),
        (SolveEvidencePhase::EdgePageRead, Some(0), 1),
        (SolveEvidencePhase::ArtifactDescriptorRead, Some(0), 0),
        (SolveEvidencePhase::ArtifactDescriptorRead, Some(0), 1),
        (SolveEvidencePhase::ImportSummaryRead, None, 1),
        (SolveEvidencePhase::ImportSummaryUtf8, None, 1),
        (SolveEvidencePhase::ImportSummaryParse, None, 1),
        (SolveEvidencePhase::RawSourceRead, Some(0), 1),
        (SolveEvidencePhase::PromotionReceiptRead, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshRead, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshPreflight, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshDecode, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshEncodeCompare, Some(0), 2),
        (SolveEvidencePhase::NamedGroupFaceRange, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportRead, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportUtf8, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportParse, Some(0), 1),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), 0),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), 1),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), u64::MAX),
        (SolveEvidencePhase::EntityResolution, None, 0),
        (SolveEvidencePhase::EntityResolution, None, 1),
        (SolveEvidencePhase::AssignmentDerivation, None, 0),
        (SolveEvidencePhase::AssignmentDerivation, None, 1),
        (SolveEvidencePhase::AssignmentDerivation, None, 2),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 0),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 1),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 2),
        (SolveEvidencePhase::ReceiptDerivation, None, 0),
        (SolveEvidencePhase::ReceiptDerivation, None, 1),
        (SolveEvidencePhase::ReceiptDerivation, None, u64::MAX),
        (SolveEvidencePhase::PrePublication, None, 0),
    ];
    assert_eq!(
        phases.len(),
        61,
        "every fresh-run solve-owned evidence phase, canonical-render entry/completion/tile boundary, duplicate-set entry/face/name traversal, and final pre-publication boundary is listed"
    );

    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let prefix_ledger = Ledger::open(":memory:").expect("prefix ledger");
    import_fixture(&prefix_ledger, &spec, bytes.clone());
    let _ = run_to_one_stage_prefix(&prefix_ledger, &decoded);
    let one_stage_prefix = solve_publication_counts(&prefix_ledger);
    for (phase, source_index, after_units) in phases {
        let ledger = Ledger::open(":memory:").expect("ledger");
        import_fixture(&ledger, &spec, bytes.clone());
        let before = solve_publication_counts(&ledger);
        let gate = CancelGate::new_clock_free();
        let plan = SolveCancellationPlan::new(phase, source_index, after_units);
        let mut clock = benign_clock();
        let mut progress = Vec::new();

        let refusal = run_solve_with_cancellation_plan(
            &ledger,
            &gate,
            &mut clock,
            &decoded,
            &mut progress,
            &plan,
        )
        .expect_err("planned in-stage cancellation refuses");
        assert!(
            plan.fired(),
            "phase {phase:?} reached its planned checkpoint"
        );
        assert_eq!(refusal.code, "cli-solve-cancelled", "phase {phase:?}");
        let has_durable_prefix =
            phase == SolveEvidencePhase::AssignmentDerivation && source_index.is_none();
        if has_durable_prefix {
            assert!(
                refusal.fix.contains("--resume"),
                "post-import phase {phase:?} recommends its durable prefix"
            );
        } else {
            assert!(
                refusal.fix.contains("retry the fresh solve"),
                "zero-prefix phase {phase:?} gives actionable fresh-retry guidance"
            );
            assert!(
                !refusal.fix.contains("--resume"),
                "zero-prefix phase {phase:?} does not recommend an unavailable checkpoint"
            );
            assert!(
                !refusal.what.contains("durable"),
                "zero-prefix phase {phase:?} does not claim a durable checkpoint"
            );
        }
        assert!(refusal.recorded_op.is_none(), "phase {phase:?}");
        assert_eq!(
            progress.len(),
            usize::from(has_durable_prefix),
            "phase {phase:?} reports exactly its durable prefix"
        );
        assert_eq!(
            solve_publication_counts(&ledger),
            if has_durable_prefix {
                one_stage_prefix
            } else {
                before
            },
            "phase {phase:?} changes no publication beyond its durable prefix"
        );
    }
}

#[test]
fn g4_chunked_raw_evidence_stops_before_a_later_storage_row_and_publishes_nothing() {
    let bytes = storage_row_spanning_tetra_stl();
    assert!(bytes.len() > STORAGE_CHUNK_LEN);
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    let imported = import_fixture(&ledger, &spec, bytes);
    let raw = imported.artifacts[0].raw_source;
    let raw_info = ledger
        .artifact_info(&raw)
        .expect("raw info query")
        .expect("raw artifact");
    assert!(raw_info.chunk_count > 0, "fixture must span storage rows");
    let before = solve_publication_counts(&ledger);

    let gate = CancelGate::new_clock_free();
    let plan = SolveCancellationPlan::new(SolveEvidencePhase::RawSourceRead, Some(0), 65_536);
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve_with_cancellation_plan(
        &ledger,
        &gate,
        &mut clock,
        &decoded,
        &mut progress,
        &plan,
    )
    .expect_err("planned chunked raw-evidence cancellation refuses");
    assert!(plan.fired());
    assert_eq!(refusal.code, "cli-solve-cancelled");
    assert!(refusal.recorded_op.is_none());
    assert!(progress.is_empty());
    assert_eq!(solve_publication_counts(&ledger), before);
}

#[test]
fn g4_exact_cap_promotion_receipt_stops_at_the_first_controlled_tile() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source, &spec, bytes);
    let entry = &imported.artifacts[0];
    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let valid_summary = String::from_utf8(
        source
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");

    let exact_cap_receipt = vec![b'r'; STORAGE_CHUNK_LEN];
    let exact_cap_hash = hash_bytes(&exact_cap_receipt);
    let exact_cap_summary = valid_summary.replacen(
        &entry.promotion_receipt.to_hex(),
        &exact_cap_hash.to_hex(),
        1,
    );
    assert_ne!(
        exact_cap_summary, valid_summary,
        "summary must bind the exact-cap receipt"
    );

    let ledger = Ledger::open(":memory:").expect("target ledger");
    reproduce_import_candidate_with_overrides(
        &source,
        &imported,
        &ledger,
        &source_row.ir,
        &exact_cap_summary,
        None,
        Some(&exact_cap_receipt),
    );
    let receipt_info = ledger
        .artifact_info(&exact_cap_hash)
        .expect("receipt descriptor query")
        .expect("exact-cap receipt");
    assert_eq!(receipt_info.len, STORAGE_CHUNK_LEN as u64);
    assert_eq!(
        receipt_info.chunk_count, 0,
        "an artifact at the storage bound remains one inline SQL row"
    );
    let before = solve_publication_counts(&ledger);

    let gate = CancelGate::new_clock_free();
    let plan = SolveCancellationPlan::new(
        SolveEvidencePhase::PromotionReceiptRead,
        Some(0),
        CONTROLLED_ARTIFACT_TILE_LEN as u64,
    );
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve_with_cancellation_plan(
        &ledger,
        &gate,
        &mut clock,
        &decoded,
        &mut progress,
        &plan,
    )
    .expect_err("planned exact-cap receipt cancellation refuses");
    assert!(plan.fired());
    assert_eq!(refusal.code, "cli-solve-cancelled");
    assert!(refusal.recorded_op.is_none());
    assert!(progress.is_empty());
    assert_eq!(solve_publication_counts(&ledger), before);
}

#[test]
fn g4_cancel_between_stages_leaves_a_durable_prefix_that_resumes_identically() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);

    // Reference: an uninterrupted run in its own ledger.
    let reference = Ledger::open(":memory:").expect("ledger");
    import_fixture(&reference, &spec, bytes.clone());
    let (reference_refusal, _) = run_to_gap(&reference, &decoded);
    let reference_run = reference_refusal.run.clone().expect("run id");

    // Interrupted: the clock's second tick (end of stage 0) requests the
    // caller-owned gate, so stage 1 observes cancellation at its boundary.
    let interrupted = Ledger::open(":memory:").expect("ledger");
    import_fixture(&interrupted, &spec, bytes);
    let gate = CancelGate::new_clock_free();
    let mut calls = 0u64;
    let gate_ref = &gate;
    let mut cancelling_clock = move || {
        calls += 1;
        if calls == 2 {
            gate_ref.request();
        }
        #[allow(clippy::cast_precision_loss)]
        {
            calls as f64 * 0.001
        }
    };
    let mut progress = Vec::new();
    let cancelled = run_solve(
        &interrupted,
        &gate,
        &mut cancelling_clock,
        &decoded,
        &mut progress,
    )
    .expect_err("cancellation refuses");
    assert_eq!(cancelled.code, "cli-solve-cancelled");
    assert!(cancelled.what.contains("1 completed stage"));
    let run = cancelled.run.clone().expect("run id");
    assert_eq!(run, reference_run, "same project, same run identity");

    // Resume with a fresh gate: the remaining real stage executes, then the
    // same gap refuses with the same owner.
    let fresh_gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut resume_progress = Vec::new();
    let resumed = resume_solve(
        &interrupted,
        &fresh_gate,
        &mut clock,
        &run,
        &mut resume_progress,
    )
    .expect_err("resume still refuses at the gap");
    assert_eq!(resumed.code, "cli-solve-stage-gap");
    assert_eq!(resumed.stage, Some("material-resolve"));
    assert_eq!(resumed.dependency, Some("frankensim-hp7tb"));

    // The interrupted-then-resumed evidence equals the uninterrupted run's:
    // identical stage receipt artifact hashes, in order.
    let reference_receipts = stage_receipt_hashes(&reference, &reference_run);
    let resumed_receipts = stage_receipt_hashes(&interrupted, &run);
    assert_eq!(reference_receipts.len(), 2);
    assert_eq!(
        reference_receipts, resumed_receipts,
        "stage evidence is bit-identical across interruption"
    );
}

#[test]
fn g4_resume_reattestation_cancellation_is_zero_publication_and_retryable() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let (gap, _) = run_to_gap(&ledger, &decoded);
    let run = gap.run.expect("run id");
    let before = solve_publication_counts(&ledger);

    let phases = [
        (SolveEvidencePhase::VisibleOpPage, Some(0), 0),
        (SolveEvidencePhase::VisibleOpPage, Some(0), 1),
        (SolveEvidencePhase::VisibleOpPage, Some(0), u64::MAX),
        (SolveEvidencePhase::CandidateOpRowRead, Some(0), 0),
        (SolveEvidencePhase::CandidateOpRowRead, Some(0), u64::MAX),
        (SolveEvidencePhase::CandidateOpTextConversion, Some(0), 0),
        (SolveEvidencePhase::CandidateOpTextConversion, Some(0), 1),
        (SolveEvidencePhase::EdgePageRead, Some(0), 0),
        (SolveEvidencePhase::EdgePageRead, Some(0), 1),
        (SolveEvidencePhase::ArtifactDescriptorRead, Some(0), 0),
        (SolveEvidencePhase::ArtifactDescriptorRead, Some(0), 1),
        (SolveEvidencePhase::EdgeSealRead, Some(0), 0),
        (SolveEvidencePhase::EdgeSealRead, Some(0), 1),
        (SolveEvidencePhase::ResumeStateRead, None, 1),
        (SolveEvidencePhase::ResumeStateDecode, None, 0),
        (SolveEvidencePhase::ResumeStateDecode, None, 1),
        (SolveEvidencePhase::ResumeProjectRead, None, 1),
        (SolveEvidencePhase::ResumeProjectUtf8, None, 1),
        (SolveEvidencePhase::ResumeProjectParse, None, 0),
        (SolveEvidencePhase::ResumeProjectParse, None, 1),
        (SolveEvidencePhase::ProjectValidation, None, 0),
        (SolveEvidencePhase::ProjectValidation, None, 1),
        (SolveEvidencePhase::ResumeProjectCanonicalCompare, None, 1),
        (SolveEvidencePhase::ProjectIdentityDerive, None, 0),
        (SolveEvidencePhase::ProjectIdentityDerive, None, 1),
        (SolveEvidencePhase::ProjectIdentityDerive, Some(0), 0),
        (SolveEvidencePhase::ProjectIdentityDerive, Some(0), 1),
        (SolveEvidencePhase::FiveExplicitsRender, None, 0),
        (SolveEvidencePhase::FiveExplicitsRender, None, 1),
        (SolveEvidencePhase::CanonicalProjectRender, None, 0),
        (SolveEvidencePhase::CanonicalProjectRender, None, 1),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), 0),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), 1),
        (SolveEvidencePhase::FiveExplicitsCompare, Some(0), u64::MAX),
        (SolveEvidencePhase::OperationContentIdentity, Some(0), 0),
        (
            SolveEvidencePhase::OperationContentIdentity,
            Some(0),
            u64::MAX,
        ),
        (SolveEvidencePhase::ResumeStageReceiptRead, None, 1),
        (SolveEvidencePhase::ResumeStageReceiptUtf8, None, 1),
        (SolveEvidencePhase::ResumeStageReceiptParse, None, 1),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(0),
            0,
        ),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(0),
            1,
        ),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(0),
            2,
        ),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(1),
            0,
        ),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(1),
            1,
        ),
        (
            SolveEvidencePhase::ResumeStageReceiptCanonicalCompare,
            Some(1),
            2,
        ),
        (SolveEvidencePhase::ImportIrParse, None, 1),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 0),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 1),
        (SolveEvidencePhase::ImportIrCanonicalCompare, None, 2),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 0),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 1),
        (SolveEvidencePhase::ImportIrDuplicateCheck, Some(0), 17),
        (SolveEvidencePhase::ImportSummaryRead, None, 1),
        (SolveEvidencePhase::ImportSummaryUtf8, None, 1),
        (SolveEvidencePhase::ImportSummaryParse, None, 1),
        (SolveEvidencePhase::RawSourceRead, Some(0), 1),
        (SolveEvidencePhase::PromotionReceiptRead, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshRead, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshPreflight, Some(0), 1),
        (SolveEvidencePhase::PromotedMeshDecode, Some(0), 2),
        (SolveEvidencePhase::PromotedMeshEncodeCompare, Some(0), 2),
        (SolveEvidencePhase::NamedGroupFaceRange, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportRead, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportUtf8, Some(0), 1),
        (SolveEvidencePhase::AssignmentReportParse, Some(0), 1),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), 0),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), 1),
        (SolveEvidencePhase::EdgeSetCompare, Some(0), u64::MAX),
        (SolveEvidencePhase::EntityResolution, None, 0),
        (SolveEvidencePhase::EntityResolution, None, 1),
        (SolveEvidencePhase::AssignmentDerivation, None, 0),
        (SolveEvidencePhase::AssignmentDerivation, None, 1),
        (SolveEvidencePhase::AssignmentDerivation, None, 2),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 0),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 1),
        (SolveEvidencePhase::AssignmentDerivation, Some(0), 2),
        (SolveEvidencePhase::ReceiptDerivation, None, 0),
        (SolveEvidencePhase::ReceiptDerivation, None, 1),
        (SolveEvidencePhase::ReceiptDerivation, None, u64::MAX),
    ];
    assert_eq!(
        phases.len(),
        79,
        "resume re-attestation lists every reachable read, UTF-8, opaque parse entry/completion, canonical-render entry/completion/tile, duplicate-set traversal, receipt entry/completion/tile, and canonical PLY phase"
    );
    for (phase, source_index, after_units) in phases {
        let gate = CancelGate::new_clock_free();
        let plan = SolveCancellationPlan::new(phase, source_index, after_units);
        let mut clock = benign_clock();
        let mut progress = Vec::new();
        let refusal = resume_solve_with_cancellation_plan(
            &ledger,
            &gate,
            &mut clock,
            &run,
            &mut progress,
            &plan,
        )
        .expect_err("planned resume re-attestation cancellation refuses");
        assert!(
            plan.fired(),
            "phase {phase:?} reached its planned checkpoint"
        );
        assert_eq!(refusal.code, "cli-solve-cancelled", "phase {phase:?}");
        assert!(refusal.recorded_op.is_none(), "phase {phase:?}");
        assert!(
            progress.is_empty(),
            "phase {phase:?} emitted no success line"
        );
        assert_eq!(
            solve_publication_counts(&ledger),
            before,
            "phase {phase:?} changed no solve publication table"
        );
    }

    let fresh_gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let retry = resume_solve(&ledger, &fresh_gate, &mut clock, &run, &mut progress)
        .expect_err("unchanged durable prefix still reaches the known gap");
    assert_eq!(retry.code, "cli-solve-stage-gap");
    assert_eq!(retry.stage, Some("material-resolve"));
    assert_eq!(
        stage_receipt_hashes(&ledger, &run).len(),
        2,
        "the normal retry preserves the complete durable stage prefix"
    );
}

#[test]
fn g4_resume_prepublication_cancellation_preserves_prefix_and_replay_identity() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);

    let reference = Ledger::open(":memory:").expect("reference ledger");
    import_fixture(&reference, &spec, bytes.clone());
    let (reference_gap, _) = run_to_gap(&reference, &decoded);
    let reference_run = reference_gap.run.expect("reference run id");
    let reference_receipts = stage_receipt_hashes(&reference, &reference_run);

    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let (prefix_refusal, _) = run_to_one_stage_prefix(&ledger, &decoded);
    let run = prefix_refusal.run.expect("run id");
    assert_eq!(run, reference_run);
    assert_eq!(stage_receipt_hashes(&ledger, &run).len(), 1);
    let before = solve_publication_counts(&ledger);

    let gate = CancelGate::new_clock_free();
    let plan = SolveCancellationPlan::new(SolveEvidencePhase::PrePublication, None, 1);
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal =
        resume_solve_with_cancellation_plan(&ledger, &gate, &mut clock, &run, &mut progress, &plan)
            .expect_err("assign-stage pre-publication cancellation refuses");
    assert!(plan.fired());
    assert_eq!(refusal.code, "cli-solve-cancelled");
    assert!(refusal.recorded_op.is_none());
    assert!(progress.is_empty());
    assert_eq!(solve_publication_counts(&ledger), before);
    assert_eq!(stage_receipt_hashes(&ledger, &run).len(), 1);

    let fresh_gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut retry_progress = Vec::new();
    let retry = resume_solve(&ledger, &fresh_gate, &mut clock, &run, &mut retry_progress)
        .expect_err("retry completes assign then reaches the known gap");
    assert_eq!(retry.code, "cli-solve-stage-gap");
    assert_eq!(
        stage_receipt_hashes(&ledger, &run),
        reference_receipts,
        "retry after pre-publication cancellation reproduces exact stage receipts"
    );
}

/// Collect stage receipt artifact hashes for a run, in stage order.
fn stage_receipt_hashes(ledger: &Ledger, run_hex: &str) -> Vec<String> {
    let run = SolveRunId::parse_hex(run_hex).expect("hex");
    let mut ids = ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .expect("ops");
    ids.sort_unstable();
    let mut receipts = Vec::new();
    for id in ids {
        let Some(row) = ledger.op(id).expect("op row") else {
            continue;
        };
        if row.session.as_deref() != Some(run.as_bytes().as_slice()) {
            continue;
        }
        if row.outcome.as_deref() != Some("ok") {
            continue;
        }
        let edges = ledger.op_artifact_edges_bounded(id, 64).expect("edges");
        for edge in &edges.edges {
            if edge.role != fs_ledger::EdgeRole::Out {
                continue;
            }
            let info = ledger
                .artifact_info(&edge.artifact)
                .expect("info")
                .expect("artifact");
            if info.kind == "solve-stage-receipt" {
                receipts.push(edge.artifact.to_hex());
            }
        }
    }
    receipts
}

#[test]
fn g5_independent_fresh_runs_retain_identical_stage_evidence() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);

    let mut all_receipts = Vec::new();
    for _ in 0..2 {
        let ledger = Ledger::open(":memory:").expect("ledger");
        import_fixture(&ledger, &spec, bytes.clone());
        let (refusal, _) = run_to_gap(&ledger, &decoded);
        let run = refusal.run.clone().expect("run id");
        all_receipts.push(stage_receipt_hashes(&ledger, &run));
    }
    assert_eq!(all_receipts[0].len(), 2);
    assert_eq!(
        all_receipts[0], all_receipts[1],
        "replay reproduces identical content identities"
    );
}

#[test]
fn g0_budget_enforcement_stops_the_run_with_an_honest_resumable_partial() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);

    // 50 wall seconds per clock call: stage 0 costs 50 s of the 60 s grant
    // (Ok, but crosses the 0.5 warning fraction); stage 1 pushes the total
    // to 100 s, which exceeds the 1.2x hard factor and pauses.
    let mut ticks = 0u64;
    let mut expensive_clock = move || {
        ticks += 1;
        #[allow(clippy::cast_precision_loss)]
        {
            ticks as f64 * 50.0
        }
    };
    let mut progress = Vec::new();
    let gate = CancelGate::new_clock_free();
    let outcome = run_solve(
        &ledger,
        &gate,
        &mut expensive_clock,
        &decoded,
        &mut progress,
    )
    .expect("budget exhaustion is an honest outcome, not a refusal");
    let SolveRunStatus::BudgetExceeded {
        resource,
        used,
        granted,
    } = outcome.status
    else {
        panic!("expected budget-exceeded, got {:?}", outcome.status);
    };
    assert_eq!(resource, "core-seconds");
    assert!(used > granted);
    assert_eq!(
        outcome.stages.len(),
        2,
        "both completed stages are reported"
    );
    assert!(
        outcome.run_receipt.is_some(),
        "the partial has a run receipt"
    );
    assert!(
        progress
            .iter()
            .any(|line| line.contains("\"status\":\"warning\"")),
        "the declared warning fraction fired before enforcement: {progress:?}"
    );

    // Resuming without more budget refuses: the recorded consumption
    // already exhausts the grant, and raising budgets changes the project
    // identity (a fresh run).
    let mut clock = benign_clock();
    let mut resume_progress = Vec::new();
    let refusal = resume_solve(
        &ledger,
        &gate,
        &mut clock,
        &outcome.run,
        &mut resume_progress,
    )
    .expect_err("resume without budget refuses");
    assert_eq!(refusal.code, "cli-solve-resume-budget");
}

#[test]
fn g3_resume_of_an_unknown_run_refuses() {
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = resume_solve(&ledger, &gate, &mut clock, &"ab".repeat(32), &mut progress)
        .expect_err("unknown run refuses");
    assert_eq!(refusal.code, "cli-solve-unknown-run");
}

#[test]
fn g3_publicly_sealed_state_without_authoritative_lineage_refuses_before_publication() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    let imported = import_fixture(&ledger, &spec, bytes);
    let run = SolveRunId::derive(&decoded);

    let project_source = ledger
        .put_artifact("solve-project-source", decoded.canonical.as_bytes(), None)
        .expect("project source");
    let forged_ir = format!(
        "{{\"schema\":\"frankensim.cli.solve-stage.v1\",\"stage\":\"assign\",\"ordinal\":1,\"run\":\"{}\",\"project\":\"{}\",\"driver_version\":{SOLVE_DRIVER_VERSION}}}",
        run.to_hex(),
        decoded.hash().to_hex(),
    );
    let forged_op = ledger
        .begin_op(
            Some(run.as_bytes()),
            &forged_ir,
            &FiveExplicits {
                seed: b"forged",
                versions: r#"{"forged":true}"#,
                budget: r#"{"forged":true}"#,
                capability: r#"{"forged":true}"#,
            },
            2,
        )
        .expect("forged op");
    ledger
        .link(forged_op, &project_source.hash, EdgeRole::In)
        .expect("forged project edge");

    // Every construction below is intentionally available through the public
    // API. Codec validity and a caller-selected session/kind must not mint
    // authority to skip both implemented stages.
    let forged_state = SolveDriverState {
        run: *run.as_bytes(),
        project: *decoded.hash().as_bytes(),
        consumed_core_s: 0.0,
        consumed_wall_s: 0.0,
        completed: vec![
            CompletedStage {
                ordinal: 0,
                op_id: forged_op,
                receipt: imported.summary_artifact,
            },
            CompletedStage {
                ordinal: 1,
                op_id: forged_op,
                receipt: imported.summary_artifact,
            },
        ],
    };
    let run_prefix: [u8; 8] = run.as_bytes()[..8].try_into().expect("run prefix");
    let envelope = LegacySnapshotV1Adapter::<SolveDriverState>::seal(
        &forged_state,
        u64::from_le_bytes(run_prefix),
    );
    let forged_checkpoint = ledger
        .put_artifact("solve-stage-state", &envelope, None)
        .expect("forged state artifact");
    ledger
        .link(forged_op, &forged_checkpoint.hash, EdgeRole::Out)
        .expect("forged state edge");
    ledger
        .finish_op(forged_op, OpOutcome::Ok, None, 3)
        .expect("finish forged op");

    let ops_before_resume = ledger.table_count("ops").expect("op count");
    let artifacts_before_resume = ledger.table_count("artifacts").expect("artifact count");
    let seals_before_resume = ledger
        .table_count("op_artifact_edge_seals")
        .expect("edge-seal count");
    let imports_before_resume = ledger.table_count("imports").expect("imports count");
    let identities_before_resume = ledger
        .table_count("op_content_identities")
        .expect("op-content-identity count");
    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = resume_solve(&ledger, &gate, &mut clock, &run.to_hex(), &mut progress)
        .expect_err("caller-mintable state has no resume authority");
    assert_eq!(refusal.code, "cli-solve-resume-identity");
    assert!(progress.is_empty(), "no forged stage may execute");
    assert_eq!(
        ledger.table_count("ops").expect("op count"),
        ops_before_resume,
        "identity refusal occurs before any stage/refusal publication"
    );
    assert_eq!(
        ledger.table_count("artifacts").expect("artifact count"),
        artifacts_before_resume,
        "identity refusal publishes no artifacts"
    );
    assert_eq!(
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("edge-seal count"),
        seals_before_resume,
        "identity refusal publishes no lineage seal"
    );
    assert_eq!(
        ledger.table_count("imports").expect("imports count"),
        imports_before_resume,
        "identity refusal mutates no import extension row"
    );
    assert_eq!(
        ledger
            .table_count("op_content_identities")
            .expect("op-content-identity count"),
        identities_before_resume,
        "identity refusal publishes no operation identity"
    );
}

#[test]
fn g3_forged_higher_checkpoint_with_substituted_predecessor_refuses_before_publication() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let run = SolveRunId::derive(&decoded);

    // Produce the canonical assign row, receipt, and checkpoint in an
    // independent reference ledger.
    let reference = Ledger::open(":memory:").expect("reference ledger");
    import_fixture(&reference, &spec, bytes.clone());
    let _ = run_to_gap(&reference, &decoded);
    let reference_stage = reference
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .expect("reference ops")
        .into_iter()
        .filter_map(|id| reference.op(id).expect("reference op row"))
        .find(|row| {
            row.session.as_deref() == Some(run.as_bytes().as_slice())
                && row.outcome.as_deref() == Some("ok")
                && row.ir.contains("\"stage\":\"assign\"")
        })
        .expect("reference assign stage");
    let reference_edges = reference
        .op_artifact_edges_bounded(reference_stage.id, 64)
        .expect("reference assign edges")
        .edges;
    let output_with_kind = |kind: &str| {
        reference_edges
            .iter()
            .find(|edge| {
                edge.role == EdgeRole::Out
                    && reference
                        .artifact_info(&edge.artifact)
                        .expect("reference artifact info")
                        .is_some_and(|info| info.kind == kind)
            })
            .map(|edge| edge.artifact)
            .unwrap_or_else(|| panic!("reference assign output kind `{kind}`"))
    };
    let reference_receipt = output_with_kind("solve-stage-receipt");
    let reference_checkpoint = output_with_kind("solve-stage-state");

    // Build the same run only through its valid first stage.
    let ledger = Ledger::open(":memory:").expect("target ledger");
    import_fixture(&ledger, &spec, bytes);
    let gate = CancelGate::new_clock_free();
    let mut calls = 0u64;
    let gate_ref = &gate;
    let mut cancelling_clock = move || {
        calls += 1;
        if calls == 2 {
            gate_ref.request();
        }
        #[allow(clippy::cast_precision_loss)]
        {
            calls as f64 * 0.001
        }
    };
    let mut progress = Vec::new();
    let cancelled = run_solve(
        &ledger,
        &gate,
        &mut cancelling_clock,
        &decoded,
        &mut progress,
    )
    .expect_err("target stops after import verification");
    assert_eq!(cancelled.code, "cli-solve-cancelled");

    // Copy the canonical higher receipt/checkpoint bytes and canonical row.
    // The one deliberate delta is the In edge: the forged operation consumes
    // its own higher checkpoint instead of the direct stage-0 predecessor.
    let copy_output = |hash, kind| {
        let bytes = reference
            .get_artifact(&hash)
            .expect("read reference output")
            .expect("reference output exists");
        ledger
            .put_artifact(kind, &bytes, None)
            .expect("copy canonical output")
    };
    let receipt = copy_output(reference_receipt, "solve-stage-receipt");
    let checkpoint = copy_output(reference_checkpoint, "solve-stage-state");
    let forged_op = ledger
        .begin_op(
            reference_stage.session.as_deref(),
            &reference_stage.ir,
            &FiveExplicits {
                seed: &reference_stage.seed,
                versions: &reference_stage.versions,
                budget: &reference_stage.budget,
                capability: &reference_stage.capability,
            },
            reference_stage.t_start,
        )
        .expect("otherwise-canonical forged stage");
    assert_eq!(
        forged_op, reference_stage.id,
        "matching histories give the forged row its canonical operation id"
    );
    ledger
        .link(forged_op, &checkpoint.hash, EdgeRole::In)
        .expect("substituted predecessor");
    ledger
        .link(forged_op, &receipt.hash, EdgeRole::Out)
        .expect("canonical receipt edge");
    ledger
        .link(forged_op, &checkpoint.hash, EdgeRole::Out)
        .expect("canonical checkpoint edge");
    ledger
        .seal_op_artifact_edges(forged_op, 3)
        .expect("exact edge-count seal");
    ledger
        .finish_op(
            forged_op,
            OpOutcome::Ok,
            None,
            reference_stage.t_end.expect("reference finish clock"),
        )
        .expect("finish otherwise-canonical forged stage");

    let counts_before = [
        ledger.table_count("ops").expect("op count"),
        ledger.table_count("artifacts").expect("artifact count"),
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("edge-seal count"),
        ledger.table_count("imports").expect("imports count"),
        ledger
            .table_count("op_content_identities")
            .expect("identity count"),
    ];
    let fresh_gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut resume_progress = Vec::new();
    let refusal = resume_solve(
        &ledger,
        &fresh_gate,
        &mut clock,
        &run.to_hex(),
        &mut resume_progress,
    )
    .expect_err("substituted predecessor has no resume authority");
    assert_eq!(refusal.code, "cli-solve-resume-identity");
    assert!(
        resume_progress.is_empty(),
        "resume executes no forged stage"
    );
    assert_eq!(
        [
            ledger.table_count("ops").expect("op count"),
            ledger.table_count("artifacts").expect("artifact count"),
            ledger
                .table_count("op_artifact_edge_seals")
                .expect("edge-seal count"),
            ledger.table_count("imports").expect("imports count"),
            ledger
                .table_count("op_content_identities")
                .expect("identity count"),
        ],
        counts_before,
        "identity refusal is zero-publication across rows, artifacts, seals, and extensions"
    );
}

#[test]
fn g3_competing_valid_longest_checkpoints_refuse_as_ambiguous() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);

    let (first, _) = run_to_gap(&ledger, &decoded);
    let (second, _) = run_to_gap(&ledger, &decoded);
    let run = first.run.expect("first run id");
    assert_eq!(second.run.as_deref(), Some(run.as_str()));

    let counts_before = [
        ledger.table_count("ops").expect("op count"),
        ledger.table_count("artifacts").expect("artifact count"),
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("edge-seal count"),
        ledger
            .table_count("op_content_identities")
            .expect("identity count"),
    ];
    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = resume_solve(&ledger, &gate, &mut clock, &run, &mut progress)
        .expect_err("equally long valid histories are ambiguous");
    assert_eq!(refusal.code, "cli-solve-resume-identity");
    assert!(refusal.what.contains("competing independently valid"));
    assert!(progress.is_empty());
    assert_eq!(
        [
            ledger.table_count("ops").expect("op count"),
            ledger.table_count("artifacts").expect("artifact count"),
            ledger
                .table_count("op_artifact_edge_seals")
                .expect("edge-seal count"),
            ledger
                .table_count("op_content_identities")
                .expect("identity count"),
        ],
        counts_before,
        "ambiguity refusal publishes nothing"
    );
}

#[test]
fn g3_import_summary_requires_the_exact_versioned_schema_not_a_hash_substring() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source_ledger = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source_ledger, &spec, bytes);
    let source_row = source_ledger
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");

    // Reproduce the accepted import writer's exact row, Five Explicits,
    // execution context, typed artifacts, and edge roles in an isolated
    // ledger. The summary schema string below is the sole admission delta.
    let ledger = Ledger::open(":memory:").expect("ledger");
    let valid_summary = String::from_utf8(
        source_ledger
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");
    let fake_summary = valid_summary.replacen(
        "frankensim.cli.geometry-import-receipt.v1",
        "alternate.geometry-import-receipt.v1",
        1,
    );
    reproduce_import_candidate(
        &source_ledger,
        &imported,
        &ledger,
        &source_row.ir,
        &fake_summary,
    );

    let (refusal, progress) =
        run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
    assert_eq!(refusal.stage, Some("import-verify"));
    assert!(
        refusal.what.contains("no completed geometry import"),
        "{refusal:?}"
    );
    assert!(
        progress.is_empty(),
        "the decoy summary must execute no stage"
    );
}

#[test]
fn g3_import_ir_requires_exact_writer_limits_and_policy_shapes() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source, &spec, bytes);
    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let valid_summary = String::from_utf8(
        source
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");

    for (label, ir) in [
        (
            "limits",
            replace_json_object_with_empty(&source_row.ir, ",\"limits\":"),
        ),
        (
            "policy",
            replace_json_object_with_empty(&source_row.ir, ",\"policy\":"),
        ),
    ] {
        let ledger = Ledger::open(":memory:").expect("isolated decoy ledger");
        reproduce_import_candidate(&source, &imported, &ledger, &ir, &valid_summary);
        let (refusal, progress) =
            run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
        assert!(
            refusal.what.contains("no completed geometry import"),
            "{label} decoy: {refusal:?}"
        );
        assert!(
            progress.is_empty(),
            "{label} decoy must execute no solve stage"
        );
    }
}

#[test]
fn g3_import_ir_duplicate_checks_preserve_order_and_bound_labels() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source, &spec, bytes);
    let (valid_gap, _) = run_to_gap(&source, &decoded);
    assert_eq!(valid_gap.stage, Some("material-resolve"));

    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let valid_summary = String::from_utf8(
        source
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");
    let canonical_group = "\"named_groups\":[{\"name\":\"fixture-all-faces\",\"faces\":[3,1,2,0]}]";
    assert!(
        source_row.ir.contains(canonical_group),
        "the accepted writer fixture preserves deliberately non-sorted face order"
    );

    let duplicate_face = source_row.ir.replacen(
        canonical_group,
        "\"named_groups\":[{\"name\":\"fixture-all-faces\",\"faces\":[3,1,1,0]}]",
        1,
    );
    let duplicate_name = source_row.ir.replacen(
        canonical_group,
        "\"named_groups\":[{\"name\":\"fixture-all-faces\",\"faces\":[3,1,2,0]},{\"name\":\"fixture-all-faces\",\"faces\":[0]}]",
        1,
    );
    let oversized_name = source_row.ir.replacen(
        "\"name\":\"fixture-all-faces\"",
        &format!("\"name\":\"{}\"", "n".repeat(4097)),
        1,
    );
    for (label, ir) in [
        ("duplicate face", duplicate_face),
        ("duplicate name", duplicate_name),
        ("oversized name", oversized_name),
    ] {
        assert_ne!(ir, source_row.ir, "{label} mutation must change the IR");
        let ledger = Ledger::open(":memory:").expect("isolated decoy ledger");
        reproduce_import_candidate(&source, &imported, &ledger, &ir, &valid_summary);
        let (refusal, progress) =
            run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
        assert!(
            refusal.what.contains("no completed geometry import"),
            "{label}: {refusal:?}"
        );
        assert!(progress.is_empty(), "{label} must execute no solve stage");
    }
}

#[test]
fn g3_import_ir_byte_caps_must_admit_the_exact_retained_raw_bytes() {
    let bytes = tetra_stl();
    assert!(bytes.len() > 1, "fixture must exceed the forged caps");
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source, &spec, bytes);
    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let valid_summary = String::from_utf8(
        source
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");

    for (label, ir) in [
        (
            "per-source",
            replace_json_unsigned(&source_row.ir, ",\"max_source_bytes\":", 1),
        ),
        (
            "aggregate",
            replace_json_unsigned(&source_row.ir, ",\"max_total_source_bytes\":", 1),
        ),
    ] {
        let ledger = Ledger::open(":memory:").expect("isolated cap ledger");
        reproduce_import_candidate(&source, &imported, &ledger, &ir, &valid_summary);
        let (refusal, progress) =
            run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
        assert!(
            refusal.what.contains("no completed geometry import"),
            "{label} cap decoy: {refusal:?}"
        );
        assert!(
            progress.is_empty(),
            "{label} cap decoy must execute no solve stage"
        );
    }
}

#[test]
fn g3_promoted_ply_with_unknown_huge_element_refuses_before_stage_publication() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let source = Ledger::open(":memory:").expect("source ledger");
    let imported = import_fixture(&source, &spec, bytes);
    assert_eq!(imported.artifacts.len(), 1, "fixture has one geometry row");
    let source_row = source
        .op(imported.op_id)
        .expect("source import op query")
        .expect("source import op");
    let entry = &imported.artifacts[0];
    let canonical_mesh = source
        .get_artifact(&entry.promoted_mesh)
        .expect("read canonical mesh")
        .expect("canonical mesh exists");
    let vertex_line = canonical_mesh
        .windows(b"element vertex ".len())
        .position(|window| window == b"element vertex ")
        .expect("writer vertex line");
    let mut hostile_mesh = Vec::with_capacity(canonical_mesh.len() + 32);
    hostile_mesh.extend_from_slice(&canonical_mesh[..vertex_line]);
    hostile_mesh.extend_from_slice(b"element junk 100000000\n");
    hostile_mesh.extend_from_slice(&canonical_mesh[vertex_line..]);
    assert!(
        hostile_mesh.len() < 1024,
        "the malformed header remains a tiny retained artifact"
    );

    let valid_summary = String::from_utf8(
        source
            .get_artifact(&imported.summary_artifact)
            .expect("read valid summary")
            .expect("valid summary exists"),
    )
    .expect("valid summary UTF-8");
    let hostile_hash = hash_bytes(&hostile_mesh);
    let hostile_summary =
        valid_summary.replacen(&entry.promoted_mesh.to_hex(), &hostile_hash.to_hex(), 1);
    assert_ne!(
        hostile_summary, valid_summary,
        "summary binds the decoy mesh"
    );

    let ledger = Ledger::open(":memory:").expect("decoy ledger");
    reproduce_import_candidate_with_mesh(
        &source,
        &imported,
        &ledger,
        &source_row.ir,
        &hostile_summary,
        Some(&hostile_mesh),
    );
    let seals_before = ledger
        .table_count("op_artifact_edge_seals")
        .expect("edge seals");
    let (refusal, progress) =
        run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-evidence");
    assert!(
        refusal.what.contains("no completed geometry import"),
        "{refusal:?}"
    );
    assert!(
        progress.is_empty(),
        "malformed retained PLY must publish no completed-stage progress"
    );
    assert_eq!(
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("edge seals"),
        seals_before,
        "malformed retained PLY must not seal a completed stage"
    );
}

#[test]
fn g3_genuine_import_above_project_solve_envelope_refuses_without_stage_publication() {
    let bytes = tetra_stl();
    let mut spec = fixture_project(7, &bytes);
    spec.budgets.as_mut().expect("budgets").memory_bytes = 1;
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let seals_before = ledger
        .table_count("op_artifact_edge_seals")
        .expect("edge seals");

    let (refusal, progress) =
        run_to_gap_expect_code(&ledger, &decoded, "cli-solve-import-envelope");
    assert_eq!(refusal.stage, Some("import-verify"));
    assert!(
        refusal.what.contains("effective solve envelope"),
        "{refusal:?}"
    );
    assert!(
        progress.is_empty(),
        "no completed stage may publish progress"
    );
    assert_eq!(
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("edge seals"),
        seals_before,
        "envelope refusal must not seal a completed stage"
    );
}

#[test]
fn g3_unrelated_wide_success_does_not_mask_an_older_valid_import() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);

    ledger.begin().expect("begin wide fixture transaction");
    let unrelated = ledger
        .begin_op(
            None,
            "{\"schema\":\"unrelated.wide-success.v1\"}",
            &FiveExplicits {
                seed: b"unrelated",
                versions: "{}",
                budget: "{}",
                capability: "{}",
            },
            0,
        )
        .expect("wide unrelated op");
    for index in 0..1025u32 {
        let artifact = ledger
            .put_artifact("unrelated-tiny-evidence", &index.to_le_bytes(), None)
            .expect("unrelated artifact");
        ledger
            .link(unrelated, &artifact.hash, EdgeRole::Out)
            .expect("unrelated edge");
    }
    ledger
        .finish_op(unrelated, OpOutcome::Ok, None, 1)
        .expect("finish unrelated op");
    ledger.commit().expect("commit wide fixture transaction");

    let (refusal, progress) = run_to_gap(&ledger, &decoded);
    assert_eq!(refusal.code, "cli-solve-stage-gap", "{refusal:?}");
    assert_eq!(refusal.stage, Some("material-resolve"));
    assert!(
        progress.iter().any(|line| line.contains("import-verify")),
        "the older valid import must still execute the solve prefix"
    );
}

fn append_unrelated_completed_ops(ledger: &Ledger, count: usize) {
    ledger.begin().expect("begin unrelated history transaction");
    for index in 0..count {
        let started = i64::try_from(index)
            .expect("bounded fixture index fits i64")
            .saturating_mul(2)
            .saturating_add(10_000);
        let op = ledger
            .begin_op(
                None,
                "{\"schema\":\"unrelated.history.v1\"}",
                &FiveExplicits {
                    seed: b"unrelated-history",
                    versions: "{}",
                    budget: "{}",
                    capability: "{}",
                },
                started,
            )
            .expect("append unrelated history op");
        ledger
            .finish_op(op, OpOutcome::Ok, None, started.saturating_add(1))
            .expect("finish unrelated history op");
    }
    ledger.commit().expect("commit unrelated history");
}

fn exact_json_bytes(len: usize) -> String {
    assert!(len >= 2);
    format!("\"{}\"", "x".repeat(len - 2))
}

#[test]
fn g3_multi_page_unrelated_history_still_finds_the_older_valid_import() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    append_unrelated_completed_ops(&ledger, 193);

    let (refusal, progress) = run_to_gap(&ledger, &decoded);
    assert_eq!(refusal.code, "cli-solve-stage-gap", "{refusal:?}");
    assert_eq!(refusal.stage, Some("material-resolve"));
    assert!(
        progress.iter().any(|line| line.contains("import-verify")),
        "descending discovery must continue across multiple pages"
    );
}

#[test]
fn g4_cap_plus_one_history_refuses_explicitly_without_publication() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    append_unrelated_completed_ops(&ledger, MAX_SOLVE_VISIBLE_OP_IDS);
    let before = solve_publication_counts(&ledger);

    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve(&ledger, &gate, &mut clock, &decoded, &mut progress)
        .expect_err("history beyond the frozen invocation cap refuses");
    assert_eq!(refusal.code, "cli-solve-work-envelope", "{refusal:?}");
    assert!(
        refusal.what.contains("visible operation ids"),
        "{refusal:?}"
    );
    assert!(refusal.recorded_op.is_none());
    assert!(progress.is_empty());
    assert_eq!(solve_publication_counts(&ledger), before);
}

#[test]
fn g4_maximum_operation_fields_are_tiled_cancelled_and_retryable() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let session = vec![b's'; MAX_OP_FIELD_BYTES];
    let seed = vec![b'r'; MAX_OP_FIELD_BYTES];
    let ir = exact_json_bytes(MAX_OP_FIELD_BYTES);
    let versions = exact_json_bytes(MAX_OP_FIELD_BYTES);
    let budget = exact_json_bytes(MAX_OP_FIELD_BYTES);
    let capability = exact_json_bytes(MAX_OP_FIELD_BYTES);
    let unrelated = ledger
        .begin_op(
            Some(&session),
            &ir,
            &FiveExplicits {
                seed: &seed,
                versions: &versions,
                budget: &budget,
                capability: &capability,
            },
            100,
        )
        .expect("maximum-field unrelated operation");
    ledger
        .finish_op(unrelated, OpOutcome::Ok, None, 101)
        .expect("finish maximum-field unrelated operation");
    let before = solve_publication_counts(&ledger);

    let gate = CancelGate::new_clock_free();
    let plan = SolveCancellationPlan::new(
        SolveEvidencePhase::CandidateOpRowRead,
        Some(0),
        u64::try_from(CONTROLLED_OP_FIELD_TILE_LEN).expect("tile size fits u64"),
    );
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = run_solve_with_cancellation_plan(
        &ledger,
        &gate,
        &mut clock,
        &decoded,
        &mut progress,
        &plan,
    )
    .expect_err("maximum-field controlled row read stops at a tile boundary");
    assert!(plan.fired());
    assert_eq!(refusal.code, "cli-solve-cancelled");
    assert!(refusal.recorded_op.is_none());
    assert!(progress.is_empty());
    assert_eq!(solve_publication_counts(&ledger), before);

    let (retry, retry_progress) = run_to_gap(&ledger, &decoded);
    assert_eq!(retry.code, "cli-solve-stage-gap", "{retry:?}");
    assert!(
        retry_progress
            .iter()
            .any(|line| line.contains("import-verify")),
        "retry must skip the unrelated maximum-field row and find the valid import"
    );
}

#[test]
fn g3_unrelated_same_run_wide_success_does_not_mask_a_valid_checkpoint() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let (initial, _) = run_to_gap(&ledger, &decoded);
    let run_hex = initial.run.expect("run id");
    let run = SolveRunId::parse_hex(&run_hex).expect("run id parses");
    let receipts_before = stage_receipt_hashes(&ledger, &run_hex);
    assert_eq!(receipts_before.len(), 2, "fixture has a valid checkpoint");

    ledger.begin().expect("begin wide fixture transaction");
    let unrelated = ledger
        .begin_op(
            Some(run.as_bytes()),
            "{\"schema\":\"unrelated.same-run-wide-success.v1\"}",
            &FiveExplicits {
                seed: b"unrelated",
                versions: "{}",
                budget: "{}",
                capability: "{}",
            },
            100,
        )
        .expect("same-run unrelated op");
    for index in 0..1025u32 {
        let artifact = ledger
            .put_artifact("unrelated-same-run-evidence", &index.to_le_bytes(), None)
            .expect("unrelated artifact");
        ledger
            .link(unrelated, &artifact.hash, EdgeRole::Out)
            .expect("unrelated edge");
    }
    ledger
        .finish_op(unrelated, OpOutcome::Ok, None, 101)
        .expect("finish unrelated op");
    ledger.commit().expect("commit wide fixture transaction");

    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let refusal = resume_solve(&ledger, &gate, &mut clock, &run_hex, &mut progress)
        .expect_err("valid checkpoint resumes to the known stage gap");
    assert_eq!(refusal.code, "cli-solve-stage-gap", "{refusal:?}");
    assert_eq!(refusal.stage, Some("material-resolve"));
    assert!(progress.is_empty(), "no new stage executes before the gap");
    assert_eq!(
        stage_receipt_hashes(&ledger, &run_hex),
        receipts_before,
        "the unrelated operation neither obscures nor extends the valid checkpoint"
    );
}

#[test]
fn g3_resume_with_a_corrupted_retained_project_refuses_identity() {
    let bytes = tetra_stl();
    let spec = fixture_project(7, &bytes);
    let decoded = decode(&spec);
    let ledger = Ledger::open(":memory:").expect("ledger");
    import_fixture(&ledger, &spec, bytes);
    let (refusal, _) = run_to_gap(&ledger, &decoded);
    let run = refusal.run.expect("run id");

    // Corrupt the retained project source; resume must refuse rather than
    // trust the damaged pin.
    let run_id = SolveRunId::parse_hex(&run).expect("hex");
    let source_hash = retained_project_source_hash(&ledger, &run_id);
    ledger
        .corrupt_artifact_for_test(&source_hash)
        .expect("test corruption");

    let gate = CancelGate::new_clock_free();
    let mut clock = benign_clock();
    let mut progress = Vec::new();
    let resumed = resume_solve(&ledger, &gate, &mut clock, &run, &mut progress)
        .expect_err("corrupted pin refuses");
    // The ledger's own read-integrity gate fires before the driver's
    // identity re-derivation gets a chance to: an earlier, equally
    // fail-closed refusal. The driver's `cli-solve-resume-identity` arm
    // remains defense-in-depth behind it.
    assert_eq!(resumed.code, "cli-solve-ledger");
    assert!(
        resumed.what.contains("reading the retained project failed"),
        "{resumed:?}"
    );
}

fn retained_project_source_hash(ledger: &Ledger, run: &SolveRunId) -> fs_ledger::ContentHash {
    let ids = ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .expect("ops");
    for id in ids {
        let Some(row) = ledger.op(id).expect("op row") else {
            continue;
        };
        if row.session.as_deref() != Some(run.as_bytes().as_slice()) {
            continue;
        }
        let edges = ledger.op_artifact_edges_bounded(id, 64).expect("edges");
        for edge in &edges.edges {
            let info = ledger
                .artifact_info(&edge.artifact)
                .expect("info")
                .expect("artifact");
            if info.kind == "solve-project-source" {
                return edge.artifact;
            }
        }
    }
    panic!("no retained project source for the run");
}

#[test]
fn g0_stage_order_and_gap_owners_are_pinned() {
    let names: Vec<&str> = SolveStage::ALL.iter().map(|stage| stage.name()).collect();
    assert_eq!(
        names,
        [
            "import-verify",
            "assign",
            "material-resolve",
            "flow-network",
            "conduction",
            "qoi",
        ],
    );
    assert_eq!(SolveStage::ImportVerify.gap_dependency(), None);
    assert_eq!(SolveStage::Assign.gap_dependency(), None);
    assert_eq!(
        SolveStage::MaterialResolve.gap_dependency(),
        Some("frankensim-hp7tb")
    );
    assert_eq!(
        SolveStage::FlowNetwork.gap_dependency(),
        Some("frankensim-frn2i")
    );
    assert_eq!(
        SolveStage::Conduction.gap_dependency(),
        Some("frankensim-s93ej")
    );
    assert_eq!(SolveStage::Qoi.gap_dependency(), Some("frankensim-s2l9v"));
}
