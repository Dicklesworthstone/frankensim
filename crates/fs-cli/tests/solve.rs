//! G0/G3/G4/G5 evidence for the solve orchestration driver
//! (bead frankensim-extreal-program-f85xj.6.5, slice 1).
//!
//! The battery drives the library seam directly: fixture project, real
//! import into an in-memory ledger, then the staged solve engine with a
//! scripted clock and caller-owned cancellation gate.

use fs_cli::{
    GeometryImportLimits, RawGeometryLibrary, SolveRefusal, SolveRunId, SolveRunStatus, SolveStage,
    import_project_geometry, resume_solve, run_solve,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::quarantine::import_mesh;
use fs_ledger::Ledger;
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
fn import_fixture(ledger: &Ledger, spec: &ProjectSpec, bytes: Vec<u8>) {
    let artifact = &spec.geometry.as_ref().expect("geometry")[0];
    let mut raw = RawGeometryLibrary::new();
    assert!(!raw.insert_mesh(
        artifact,
        "fixtures/enclosure.stl",
        bytes,
        "m",
        0,
        Vec::new(),
    ));
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        import_project_geometry(spec, &raw, ledger, GeometryImportLimits::DEFAULT, cx)
            .expect("fixture imports")
    });
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

#[test]
fn g0_run_identity_is_deterministic_and_input_sensitive() {
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
