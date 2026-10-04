//! G0 command-contract tests for the `frankensim` CLI membrane.

use fs_cli::{
    MAX_CARD_PACK_BYTES, MAX_CARD_PACK_SOURCE_BYTES, MAX_CARD_PACKS, exit, run, validate_source,
};
use fs_project::{
    Budgets, ConsequenceClass, Cooling, DecisionGate, EntityDecl, Envelope, GeometryArtifact,
    GeometryAssignment, MeshSelector, Metadata, OutputRequest, PowerDissipation, ProjectSpec,
    RequirementDirection, RequirementSeverity, RequirementSource, RequirementSourceKind,
    SafetyFactorPolicy, Seeds, SolverSettings, ThermalLimit, UnitsDoctrine, Versions, print_sexpr,
};
use fs_qty::QtyAny;

fn args(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

/// Per-test scratch directory under the platform temp root.
///
/// The card-pack resource ceilings are the only part of the CLI contract that
/// has to reach real filesystem metadata — a size ceiling that never sees a
/// real `stat` is not a ceiling — so this is the one place the membrane tests
/// touch disk.
fn scratch(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("fs-cli-cards-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

/// Write the admissible fixture project so a solve invocation gets past the
/// project read and reaches card-pack admission.
fn written_project(dir: &std::path::Path) -> std::path::PathBuf {
    let source = print_sexpr(&valid_project()).expect("fixture renders canonically");
    let path = dir.join("reference.fsim");
    std::fs::write(&path, source).expect("project fixture writes");
    path
}

fn solve_with_pack(project: &std::path::Path, ledger: &std::path::Path, pack: &str) -> Vec<String> {
    vec![
        "solve".to_string(),
        project.to_string_lossy().into_owned(),
        ledger.to_string_lossy().into_owned(),
        "--materials".to_string(),
        pack.to_string(),
        "--json".to_string(),
    ]
}

fn valid_project() -> ProjectSpec {
    let kelvin = |value| QtyAny::new(value, fs_project::spec::dims::TEMPERATURE);
    let watts = |value| QtyAny::new(value, fs_project::spec::dims::POWER);
    ProjectSpec {
        metadata: Some(Metadata {
            name: "cli-reference".to_string(),
            created: "2026-07-22".to_string(),
            context_of_use: "CLI contract conformance".to_string(),
            intended_decision: "exercise structural project admission".to_string(),
            decision_gate: DecisionGate::ScopingEstimate,
            consequence: ConsequenceClass::Advisory,
        }),
        versions: Some(Versions {
            schema: fs_project::FSIM_VERSION,
            constellation: "00".repeat(32),
            workspace: "11".repeat(20),
        }),
        seeds: Some(Seeds { root: 7 }),
        budgets: Some(Budgets {
            solve_time: QtyAny::new(60.0, fs_project::spec::dims::TIME),
            memory_bytes: 1024 * 1024,
            accuracy_rel: 0.01,
        }),
        capabilities: Some(vec!["thermal.conduction-solve".to_string()]),
        units: Some(UnitsDoctrine {
            storage: "si-base".to_string(),
            display: "engineering".to_string(),
        }),
        geometry: Some(vec![GeometryArtifact {
            role: "plate".to_string(),
            format: "stl".to_string(),
            source_hash: 9,
            parser_version: "1".to_string(),
            surface_offset: None,
        }]),
        assignments: Some(vec![GeometryAssignment {
            artifact: "plate".to_string(),
            target: "hot".to_string(),
            length_unit: "m".to_string(),
            selector: MeshSelector::NamedGroup {
                name: "HOT".to_string(),
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
                name: "plate".to_string(),
                display: "Plate".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "plate".to_string(),
                name: "hot".to_string(),
                display: "Hot region".to_string(),
                expect_id: None,
            },
        ]),
        materials: Some(Vec::new()),
        interface_cards: Some(Vec::new()),
        perfect_contacts: None,
        power: Some(vec![PowerDissipation {
            region: "hot".to_string(),
            watts: watts(5.0),
            duty: 1.0,
        }]),
        cooling: Some(Cooling {
            fans: Vec::new(),
            vents: Vec::new(),
            leakage: watts(0.0),
            airflow_leakage: None,
            fan_system: None,
            conduction: None,
            fan_efficiency: None,
        }),
        envelope: Some(Envelope {
            ambient_lo: kelvin(293.15),
            ambient_hi: kelvin(313.15),
            pressure: QtyAny::new(101_325.0, fs_project::spec::dims::PRESSURE),
        }),
        requirements: Some(vec![ThermalLimit {
            qoi: "temperature-max".to_string(),
            class: "surface".to_string(),
            region: "hot".to_string(),
            direction: RequirementDirection::AtMost,
            limit: kelvin(353.15),
            margin: kelvin(5.0),
            source: RequirementSource {
                kind: RequirementSourceKind::UserDeclaration,
                document: "cli-test-declaration".to_string(),
                version: "1".to_string(),
                locator: "temperature-max".to_string(),
            },
            safety_factor: SafetyFactorPolicy {
                factor: 1.0,
                source: RequirementSource {
                    kind: RequirementSourceKind::UserDeclaration,
                    document: "cli-test-margin-policy".to_string(),
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
            region: None,
        }]),
    }
}

#[test]
fn g0_validate_accepts_only_a_strictly_admissible_project() {
    let source = print_sexpr(&valid_project()).expect("fixture renders canonically");
    let output = validate_source("reference.fsim", &source, false, true);
    assert_eq!(output.exit_code, exit::SUCCESS);
    assert!(output.stderr.is_empty());
    assert!(output.stdout.contains("\"status\":\"ok\""));
    assert!(output.stdout.contains("\"finding_count\":0"));
    assert!(
        output
            .stdout
            .contains("\"authority\":\"structural-project-admission\"")
    );
    assert_eq!(output.stdout.lines().count(), 1, "one JSON result record");
}

#[test]
fn g0_validate_retains_every_finding_and_fix() {
    let source = print_sexpr(&ProjectSpec::default()).expect("draft renders");
    let output = validate_source("broken.fsim", &source, false, true);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stdout.contains("\"status\":\"refused\""));
    assert!(output.stdout.contains("\"finding_count\":17"));
    assert_eq!(output.stderr.lines().count(), 17);
    assert!(output.stderr.contains("project-metadata-missing"));
    assert!(output.stderr.contains("\"fix\":"));
}

#[test]
fn g0_validate_refuses_noncanonical_bytes_without_rewriting_them() {
    let mut source = print_sexpr(&valid_project()).expect("fixture renders");
    source.push('\n');
    let output = validate_source("reference.fsim", &source, false, false);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stderr.contains("fsim-non-canonical"));
    assert!(output.stderr.contains("use the lenient parser"));
}

#[test]
fn g0_argument_grammar_and_json_flag_are_stable() {
    let help = run(args(&["--json", "help"]));
    assert_eq!(help.exit_code, exit::SUCCESS);
    assert!(help.stdout.contains("\"command\":\"help\""));
    assert!(
        help.stdout
            .contains("import <project> <source>... <ledger.db>")
    );

    let duplicate = run(args(&["validate", "x.fsim", "--json", "--json"]));
    assert_eq!(duplicate.exit_code, exit::USAGE);
    assert!(duplicate.stderr.contains("cli-duplicate-flag"));

    let extra = run(args(&["report", "run-1", "ledger.db", "extra"]));
    assert_eq!(extra.exit_code, exit::USAGE);
    assert!(extra.stderr.contains("cli-usage"));

    let unknown_flag = run(args(&["validate", "--lenient"]));
    assert_eq!(unknown_flag.exit_code, exit::USAGE);
    assert!(unknown_flag.stderr.contains("cli-usage"));

    assert!(
        help.stdout.contains("[--materials <pack>]"),
        "the published usage names the card-pack grammar"
    );

    let mixed_import_policy = run(args(&[
        "import",
        "project.fsim",
        "mesh.stl",
        "run.db",
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
        "--step-root",
        "60",
        "--target-h",
        "1",
    ]));
    assert_eq!(mixed_import_policy.exit_code, exit::USAGE);
    assert!(mixed_import_policy.stderr.contains("cli-import-usage"));

    let invalid_spacing = run(args(&[
        "import",
        "project.fsim",
        "mesh.step",
        "run.db",
        "--unit",
        "m",
        "--step-root",
        "60",
        "--target-h",
        "NaN",
    ]));
    assert_eq!(invalid_spacing.exit_code, exit::USAGE);
    assert!(invalid_spacing.stderr.contains("cli-import-argument"));
}

#[test]
fn g0_report_and_package_refuse_an_unknown_run_without_writing_anything() {
    // The export verbs read only what a completed solve retained. A ledger
    // that never saw the run must yield the solve loader's own refusal code,
    // and no report, twin, or package file may appear on disk.
    let dir = scratch("typed-stage-gaps");
    let ledger = dir.join("fixture-ledger.db");
    let _ = fs_ledger::Ledger::open(ledger.to_str().expect("UTF-8 fixture path"))
        .expect("fixture ledger opens");
    let run_id = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

    for verb in ["report", "package"] {
        let output = run(vec![
            "--json".to_string(),
            verb.to_string(),
            run_id.to_string(),
            ledger.to_string_lossy().into_owned(),
        ]);
        assert_eq!(output.exit_code, exit::REFUSED, "{verb}: {}", output.stderr);
        assert!(
            output.stderr.contains("cli-solve-unknown-run"),
            "{verb}: {}",
            output.stderr
        );
        assert!(output.stdout.contains("\"status\":\"refused\""));
        assert!(output.stdout.contains(&format!("\"command\":\"{verb}\"")));
        assert!(!output.stdout.contains("342.15"));
        assert!(!output.stdout.contains("\"merkle_root\""));
        assert!(!output.stdout.contains("\"content_hash\""));
    }
    for suffix in [".report.html", ".report.json", ".fspkg"] {
        assert!(
            !dir.join(format!("{run_id}{suffix}")).exists(),
            "a refused export must not write {suffix}"
        );
    }
}

#[test]
fn g0_solve_grammar_requires_a_ledger_operand() {
    // The pre-6.5 one-operand spellings are now off-grammar, not unavailable.
    let missing_ledger = run(args(&["solve", "project.fsim", "--json"]));
    assert_eq!(missing_ledger.exit_code, exit::USAGE);
    assert!(missing_ledger.stderr.contains("cli-usage"));

    let missing_resume_ledger = run(args(&["solve", "--resume", "run-1", "--json"]));
    assert_eq!(missing_resume_ledger.exit_code, exit::USAGE);
    assert!(missing_resume_ledger.stderr.contains("cli-usage"));

    // A well-formed solve against a missing project fails at bounded input,
    // before any ledger side effect.
    let missing_project = run(args(&["solve", "no-such.fsim", "no-such.db", "--json"]));
    assert_eq!(missing_project.exit_code, exit::INPUT);
    assert!(missing_project.stderr.contains("cli-input-read"));
}

#[test]
fn g0_solve_card_pack_flags_are_repeatable_and_pair_strictly() {
    // A dangling flag with no value is a usage refusal, not a silent drop.
    let dangling = run(args(&[
        "solve",
        "no-such.fsim",
        "no-such.db",
        "--materials",
        "--json",
    ]));
    assert_eq!(dangling.exit_code, exit::USAGE);
    assert!(dangling.stderr.contains("cli-solve-usage"));

    let unknown = run(args(&[
        "solve",
        "no-such.fsim",
        "no-such.db",
        "--cards",
        "p.fsmcdpk",
        "--json",
    ]));
    assert_eq!(unknown.exit_code, exit::USAGE);
    assert!(unknown.stderr.contains("cli-solve-usage"));

    // Repetition is legal grammar: the project reads the missing-input
    // refusal, which proves parsing accepted both pairs and got as far as
    // bounded project I/O.
    let repeated = run(args(&[
        "solve",
        "no-such.fsim",
        "no-such.db",
        "--materials",
        "a.fsmcdpk",
        "--materials",
        "b.fsmcdpk",
        "--interfaces",
        "c.fsintpk",
        "--json",
    ]));
    assert_eq!(repeated.exit_code, exit::INPUT);
    assert!(repeated.stderr.contains("cli-input-read"));

    // The resume spelling keeps its own exact arity and is not reinterpreted
    // as a project/ledger pair with trailing flags.
    let resume_with_cards = run(args(&[
        "solve",
        "--resume",
        "run-1",
        "run.db",
        "--materials",
        "a.fsmcdpk",
    ]));
    assert_eq!(resume_with_cards.exit_code, exit::USAGE);
    assert!(resume_with_cards.stderr.contains("cli-usage"));
}

#[test]
fn g0_the_invocation_card_pack_ceiling_refuses_before_any_file_is_touched() {
    // The grammar ceiling counts declared pairs, so it must fire during
    // argument parsing — before the project path is even stat'd. The missing
    // project is the positive control: at the ceiling the run gets far enough
    // to fail on it, one pair past the ceiling it never does.
    let invocation = |pairs: usize| {
        let mut argv = vec![
            "solve".to_string(),
            "no-such.fsim".to_string(),
            "no-such.db".to_string(),
        ];
        for index in 0..pairs {
            argv.push("--materials".to_string());
            argv.push(format!("pack-{index}.fsmcdpk"));
        }
        argv.push("--json".to_string());
        run(argv)
    };

    let at_cap = invocation(MAX_CARD_PACKS);
    assert_eq!(at_cap.exit_code, exit::INPUT);
    assert!(
        at_cap.stderr.contains("cli-input-read"),
        "exactly the ceiling parses and proceeds to bounded project I/O"
    );

    let past_cap = invocation(MAX_CARD_PACKS + 1);
    assert!(past_cap.stderr.contains("cli-solve-card-pack-count"));
    assert!(
        !past_cap.stderr.contains("cli-input-read"),
        "the ceiling must refuse before the project is read, not after"
    );
    // Endorsed (bead p63op): the invocation matches the documented grammar —
    // `[--materials <pack>]...` is unbounded repetition — and is refused by a
    // resource ceiling, which is exactly the case the `exit::INPUT` doc names.
    // The same code now reaches this class from every layer that can emit it.
    assert_eq!(past_cap.exit_code, exit::INPUT);
}

#[test]
fn g0_a_non_regular_card_pack_path_refuses_at_the_size_guard() {
    let dir = scratch("nonregular");
    let project = written_project(&dir);
    // A directory resolves through `stat` but is not a regular file, so it
    // can never carry a bounded pack read.
    let output = run(solve_with_pack(
        &project,
        &dir.join("run.db"),
        &dir.to_string_lossy(),
    ));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stderr.contains("cli-solve-card-pack-size"));
}

#[test]
fn g0_the_card_pack_read_ceiling_is_exactly_max_card_pack_bytes_on_disk() {
    let dir = scratch("oversized");
    let project = written_project(&dir);
    let ledger = dir.join("run.db");

    // Both files are undecodable, so the code is what discriminates: one byte
    // past the ceiling never reaches the decoder, exactly at the ceiling does.
    let past_cap = dir.join("past-cap.fsmcdpk");
    std::fs::write(&past_cap, vec![0u8; MAX_CARD_PACK_BYTES as usize + 1])
        .expect("oversized fixture writes");
    let output = run(solve_with_pack(
        &project,
        &ledger,
        &past_cap.to_string_lossy(),
    ));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stderr.contains("cli-solve-card-pack-size"));

    let at_cap = dir.join("at-cap.fsmcdpk");
    std::fs::write(&at_cap, vec![0u8; MAX_CARD_PACK_BYTES as usize])
        .expect("at-ceiling fixture writes");
    let output = run(solve_with_pack(
        &project,
        &ledger,
        &at_cap.to_string_lossy(),
    ));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(
        output.stderr.contains("cli-solve-card-pack-decode"),
        "bytes exactly at the ceiling are read in full and refused by the decoder"
    );
}

#[test]
fn g0_an_overlong_pack_path_refuses_as_unreadable_not_as_an_oversized_label() {
    // `cli-solve-card-pack-source` guards the retained diagnostic label, but
    // from the CLI the label IS the path, and no filesystem admits a
    // component this long. The read guard therefore shadows it: the source
    // ceiling is a library-boundary guard only, proven reachable in
    // `fs_cli::cards`' own unit battery rather than pretended to be covered
    // here. If the guard order ever changes, this pin is what says so.
    let dir = scratch("longpath");
    let project = written_project(&dir);
    let overlong = dir.join("x".repeat(MAX_CARD_PACK_SOURCE_BYTES + 1));
    let output = run(solve_with_pack(
        &project,
        &dir.join("run.db"),
        &overlong.to_string_lossy(),
    ));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stderr.contains("cli-solve-card-pack-read"));
    assert!(!output.stderr.contains("cli-solve-card-pack-source"));
}

#[test]
fn g0_json_diagnostics_escape_user_controlled_subjects() {
    let output = validate_source("bad\"name\n.fsim", "not a project", false, true);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stderr.contains("bad\\\"name\\n.fsim"));
    assert_eq!(output.stderr.lines().count(), 1);
}

#[test]
fn g0_validate_path_refuses_unknown_extensions_before_reading() {
    let output = run(args(&["validate", "project.yaml", "--json"]));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stderr.contains("cli-input-format"));
    assert!(output.stderr.contains(".fsim or .json"));
}

#[test]
fn g0_import_command_routes_valid_policy_to_bounded_project_io() {
    for invocation in [
        &[
            "import",
            "missing.fsim",
            "mesh.stl",
            "run.db",
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ][..],
        &[
            "import",
            "missing.fsim",
            "mesh.step",
            "run.db",
            "--unit",
            "m",
            "--step-root",
            "60",
            "--target-h",
            "1",
        ][..],
    ] {
        let output = run(args(invocation));
        assert_eq!(output.exit_code, exit::INPUT);
        assert!(output.stdout.contains("command=import"));
        assert!(output.stderr.contains("cli-input-read"));
    }
}

#[test]
fn g0_the_tracked_reference_project_validates_through_the_real_cli_verb() {
    // Every other project in this battery is built in-process. This one is
    // read off disk through the actual product verb, which is the only way
    // to prove the documented user story ("write a .fsim, validate it")
    // has a starting point that works (bead frankensim-58fbi).
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/reference-project/cooling-reference.fsim");
    let output = run(args(&["--json", "validate", &path.to_string_lossy()]));
    assert_eq!(output.exit_code, exit::SUCCESS, "stderr: {}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"ok\""));
    assert!(output.stdout.contains("\"finding_count\":0"));
    assert_eq!(output.stdout.lines().count(), 1, "one JSON result record");
}

#[test]
fn g0_the_worked_example_fixtures_stay_fresh_through_the_real_cli_verb() {
    // The worked examples (bead frankensim-extreal-program-f85xj.6.12) are
    // executed, not prose. The minimal heated-plate fixture must keep
    // validating clean; the refusal-loop fixture must keep refusing with
    // exactly the documented code; and its one-token repair must remain
    // byte-identical to the tracked reference project.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");

    let heated = root.join("examples/heated-plate/heated-plate.fsim");
    let output = run(args(&["--json", "validate", &heated.to_string_lossy()]));
    assert_eq!(output.exit_code, exit::SUCCESS, "stderr: {}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"ok\""));
    assert!(output.stdout.contains("\"finding_count\":0"));
    // Frozen current-schema canonical hashes: the real verb migrates these
    // historical fixtures before hashing. Physical fixture drift still fails
    // here; an intentional schema migration updates these pins using the
    // real verb while retaining the original fixture bytes.
    assert!(
        output.stdout.contains(
            "\"project_hash\":\"d3b7c322ee44da3f165dd07902d3c6d943ba26b4a475464f8c43486b0c32cf1b\""
        ),
        "heated-plate.fsim drifted from its frozen canonical hash"
    );

    let reference = root.join("data/reference-project/cooling-reference.fsim");
    let ref_out = run(args(&["--json", "validate", &reference.to_string_lossy()]));
    assert_eq!(ref_out.exit_code, exit::SUCCESS);
    assert!(
        ref_out.stdout.contains(
            "\"project_hash\":\"6c6eb9783a8fe6a1a0278e51d658c592f505a907b6a71d3c569165244297b863\""
        ),
        "cooling-reference.fsim drifted from its frozen canonical hash"
    );

    let broken = root.join("examples/refusal-loop/broken.fsim");
    let output = run(args(&["--json", "validate", &broken.to_string_lossy()]));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(
        output.stderr.contains("project-duty-range"),
        "stderr: {}",
        output.stderr
    );
    assert!(
        output.stderr.contains("duty must lie in 0.0..=1.0"),
        "stderr: {}",
        output.stderr
    );

    let reference_bytes = std::fs::read(root.join("data/reference-project/cooling-reference.fsim"))
        .expect("tracked reference project is readable");
    let broken_text =
        std::fs::read_to_string(&broken).expect("refusal-loop fixture is readable utf-8");
    let repaired = broken_text.replacen(":duty 2.0", ":duty 1.0", 1);
    assert_eq!(
        repaired.as_bytes(),
        reference_bytes.as_slice(),
        "the one-token repair must reproduce the tracked reference bytes"
    );
}

#[test]
fn g1_the_heatsink_fan_example_runs_every_stage_through_the_real_cli_verb() {
    // The heatsink+fan worked example (bead frankensim-extreal-program-
    // f85xj.6.12; conduction declared under rc-root-q61wp.8) is the deepest
    // walkthrough the product supports: a real finned body (one closed
    // 108-facet shell), a declared fan system, vent, airflow-leakage bypass,
    // a seeded aluminium region, and an airflow-convection boundary whose
    // coefficient is derived from the vent branch's operating point through
    // the Hausen developing-flow card. It must clear all seven solve stages
    // through the one-command `run` verb and export a report and package
    // next to the ledger, every value of which traces to a retained receipt.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");

    let validated = run(args(&[
        "--json",
        "validate",
        fsim.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        validated.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        validated.stderr
    );
    assert!(validated.stdout.contains("\"status\":\"ok\""));
    assert!(validated.stdout.contains("\"finding_count\":0"));

    let dir = scratch("heatsink-run");
    let ledger = dir.join("heatsink.db");
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );

    let output = run(args(&[
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    assert!(output.stdout.contains("\"status\":\"completed\""));
    assert!(output.stdout.contains("\"stages_completed\":7"));
    // Conduction executed (not a typed gap), and the retained verdict is the
    // honest Estimated/indeterminate one with a checker-clean package.
    assert!(
        output
            .stderr
            .contains("\"stage\":\"conduction\",\"ordinal\":4,\"status\":\"ok\""),
        "stderr: {}",
        output.stderr
    );
    // This single-rung project has no discretization estimate, so one term
    // stays NO-DATA and the verdict stays Estimated/indeterminate. The ladder
    // variant measures all eight (see scripts/ci/examples_freshness_e2e.sh).
    assert!(
        output.stdout.contains("\"verdict\":\"indeterminate\""),
        "stdout: {}",
        output.stdout
    );
    assert!(
        output.stdout.contains("\"checker\":\"pass\""),
        "stdout: {}",
        output.stdout
    );
    let run_id = output
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.get(..64))
        .expect("run result names its 64-hex run id");
    assert!(
        run_id.chars().all(|c| c.is_ascii_hexdigit()),
        "run id {run_id}"
    );
    for suffix in [".report.html", ".report.json", ".fspkg"] {
        let path = dir.join(format!("{run_id}{suffix}"));
        let bytes = std::fs::read(&path)
            .unwrap_or_else(|error| panic!("{} was not exported: {error}", path.display()));
        assert!(!bytes.is_empty(), "{} is empty", path.display());
    }
    // A repeat solve of identical inputs re-attests the retained run instead
    // of driving a second chain: before 2026-09-29 the second chain (same
    // run, different wall seconds) made every later export refuse with
    // `cli-solve-resume-identity` on two competing complete checkpoints.
    let first_receipt = output
        .stdout
        .split("\"run_receipt\":\"")
        .nth(1)
        .map(|rest| rest[..64].to_string());
    let again = run(args(&[
        "--json",
        "solve",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(again.exit_code, exit::SUCCESS, "stderr: {}", again.stderr);
    assert!(
        again.stdout.contains(&format!("\"run\":\"{run_id}\"")),
        "{}",
        again.stdout
    );
    assert!(
        again.stdout.contains("\"stages_completed\":7"),
        "{}",
        again.stdout
    );
    if let Some(receipt) = &first_receipt {
        assert!(
            again.stdout.contains(receipt.as_str()),
            "re-attest must name the retained run receipt: {}",
            again.stdout
        );
    }
    let exported = run(args(&[
        "--json",
        "report",
        run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        exported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        exported.stderr
    );
    let rerun = run(args(&[
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(rerun.exit_code, exit::SUCCESS, "stderr: {}", rerun.stderr);
    assert!(
        rerun.stdout.contains("\"checker\":\"pass\""),
        "{}",
        rerun.stdout
    );
    // The declared surface offset moves T_max the physical way: more wetted
    // area and thicker fins (outward) cool the part at fixed watts, and the
    // inward bound heats it.
    let report = std::fs::read_to_string(dir.join(format!("{run_id}.report.json"))).unwrap();
    let solved = |label: &str| -> f64 {
        let key = format!("declared surface offset {label} = ");
        let tail = &report[report
            .find(&key)
            .unwrap_or_else(|| panic!("no {label} vertex in {report}"))
            + key.len()..];
        tail[..tail.find(' ').unwrap()].parse().unwrap()
    };
    let (inward, outward) = (solved("inward"), solved("outward"));
    assert!(
        outward < inward,
        "outward {outward} K must be cooler than inward {inward} K"
    );
}

#[test]
fn g1_chip_footprint_power_enters_through_the_declared_surface() {
    // q61wp.52: the heatsink with its 3 W entering through a 20 x 20 mm die
    // contact (an fsim v8 `(surface ...)` entity) instead of volumetrically.
    // The power row names the surface, so solve must lower it to an inward
    // Neumann flux over exactly the footprint (no volumetric source) and
    // carve the footprint out of the convection row; the balance then closes
    // 3 W in through the chip against 3 W out by convection. Falsifiers: the
    // measured patch area is the declared 400 mm^2, and moving the footprint
    // (same STL, box shifted -20 mm in y onto another whole-facet patch of
    // the generator's grid) moves the hot spot with it.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-fan-chip.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink-chip.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("heatsink-chip");
    let solve_chip = |project: &std::path::Path, tag: &str| -> String {
        let ledger = dir.join(format!("{tag}.db"));
        let imported = run(args(&[
            "--json",
            "import",
            project.to_string_lossy().as_ref(),
            stl.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]));
        assert_eq!(
            imported.exit_code,
            exit::SUCCESS,
            "stderr: {}",
            imported.stderr
        );
        let solved = run(args(&[
            "--json",
            "solve",
            project.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--materials",
            pack.to_string_lossy().as_ref(),
        ]));
        assert_eq!(
            solved.exit_code,
            exit::SUCCESS,
            "stdout: {} / stderr: {}",
            solved.stdout,
            solved.stderr
        );
        assert!(
            solved.stdout.contains("\"stages_completed\":7"),
            "stdout: {}",
            solved.stdout
        );
        let run_id = solved
            .stdout
            .split("\"run\":\"")
            .nth(1)
            .and_then(|rest| rest.get(..64))
            .expect("solve result names its run id");
        let ledger = fs_ledger::Ledger::open(ledger.to_str().expect("utf-8 path")).expect("ledger");
        let receipts = stage_receipt_hashes(&ledger, run_id);
        let text = receipt_text(&ledger, &receipts[4]);
        assert!(text.contains("\"stage\":\"conduction\""), "{text}");
        text
    };
    let hottest = |text: &str| -> [f64; 3] {
        let tail = text
            .split("\"hottest_vertex_m\":[")
            .nth(1)
            .unwrap_or_else(|| panic!("no hottest vertex in {text}"));
        let mut xyz = tail[..tail.find(']').unwrap()]
            .split(',')
            .map(|v| v.parse::<f64>().unwrap());
        [
            xyz.next().unwrap(),
            xyz.next().unwrap(),
            xyz.next().unwrap(),
        ]
    };

    let text = solve_chip(&fsim, "centre");
    let field = |key: &str| number_after(&text, &format!("\"{key}\":"));
    assert_eq!(
        field("source_w"),
        0.0,
        "surface power must not become a volumetric source"
    );
    assert!(
        (field("neumann_out_w") + 3.0).abs() < 1e-9,
        "3 W must enter through the footprint"
    );
    assert!(
        (field("robin_out_w") - 3.0).abs() < 1e-6,
        "3 W must leave by convection"
    );
    assert!(field("relative_closure").abs() < 1e-6);
    // The STL carries six-decimal coordinates; the area is exact to roundoff.
    assert!((field("area_m2") / 4e-4 - 1.0).abs() < 1e-6, "{text}");
    let centre = hottest(&text);
    assert!(
        (0.02..=0.04).contains(&centre[1]) && centre[2] < 1e-6,
        "{centre:?}"
    );

    let declared = std::fs::read_to_string(&fsim).unwrap();
    let footprint = "(box :min (vec3 0.03 0.02 -1e-6) :max (vec3 0.05 0.04 1e-6)";
    assert_eq!(declared.matches(footprint).count(), 1);
    let moved = dir.join("heatsink-fan-chip-moved.fsim");
    std::fs::write(
        &moved,
        declared.replace(
            footprint,
            "(box :min (vec3 0.03 0.0 -1e-6) :max (vec3 0.05 0.02 1e-6)",
        ),
    )
    .unwrap();
    let text = solve_chip(&moved, "moved");
    assert!(
        (number_after(&text, "\"area_m2\":") / 4e-4 - 1.0).abs() < 1e-6,
        "{text}"
    );
    let shifted = hottest(&text);
    assert!(
        shifted[1] <= 0.02 && shifted[2] < 1e-6,
        "hot spot {shifted:?} did not follow the footprint"
    );
    eprintln!("chip hot spot: centre {centre:?} -> moved {shifted:?}");
    assert!(
        (shifted[0] - centre[0]).abs() < 0.01,
        "x barely moves: {centre:?} -> {shifted:?}"
    );
}

#[test]
fn g1_the_contact_pair_imports_two_sources_and_conducts_through_the_declared_joint() {
    // q61wp.49/.53: two bodies, one declared card-backed contact joint, through
    // the product verbs. `import` binds one source per geometry row in
    // declaration order. All 5 W generated in the hot body must cross the
    // joint into the fixed-temperature cold body (its outer faces carry zero
    // flux), so conservation fixes the contact heat at -5 W and the mean jump
    // at -0.5 K for R'' = 0.1 m^2 K/W on the unit joint, whatever the mesh.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let example = root.join("examples/contact-pair");
    let fsim = example.join("contact-pair.fsim");
    let cold = example.join("cold-body.stl");
    let hot = example.join("hot-body.stl");
    let dir = scratch("contact-pair");
    let import = |first: &std::path::Path, second: &std::path::Path, name: &str| {
        run(args(&[
            "--json",
            "import",
            fsim.to_string_lossy().as_ref(),
            first.to_string_lossy().as_ref(),
            second.to_string_lossy().as_ref(),
            dir.join(name).to_string_lossy().as_ref(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]))
    };
    // Sources bind by declaration order; swapped bytes fail the pinned hash.
    let swapped = import(&hot, &cold, "swapped.db");
    assert_eq!(
        swapped.exit_code,
        exit::REFUSED,
        "stdout: {}",
        swapped.stdout
    );
    assert!(
        swapped.stderr.contains("cli-import-source-hash-mismatch"),
        "{}",
        swapped.stderr
    );
    let one = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        cold.to_string_lossy().as_ref(),
        dir.join("one.db").to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert!(
        one.stderr.contains("cli-import-source-count"),
        "{}",
        one.stderr
    );

    let imported = import(&cold, &hot, "pair.db");
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );
    assert!(
        imported.stdout.contains("\"artifact_count\":2"),
        "{}",
        imported.stdout
    );
    let ledger = dir.join("pair.db");
    let solved = run(args(&[
        "--json",
        "solve",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        root.join("data/reference-project/aa6061.fsmcdpk")
            .to_string_lossy()
            .as_ref(),
        "--interfaces",
        example.join("cold-hot.fsintpk").to_string_lossy().as_ref(),
    ]));
    assert_eq!(solved.exit_code, exit::SUCCESS, "stderr: {}", solved.stderr);
    assert!(
        solved.stdout.contains("\"stages_completed\":7"),
        "{}",
        solved.stdout
    );
    let run_id = solved
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.get(..64))
        .expect("run id");
    let ledger = fs_ledger::Ledger::open(ledger.to_str().unwrap()).expect("ledger");
    let receipts = stage_receipt_hashes(&ledger, run_id);
    let text = receipt_text(&ledger, &receipts[4]);
    let field = |key: &str| number_after(&text, &format!("\"{key}\":"));
    assert!(text.contains("\"interface\":\"cold-hot-joint\""), "{text}");
    assert!((field("source_w") - 5.0).abs() < 1e-11, "{text}");
    assert!((field("heat_rate_a_to_b_w") + 5.0).abs() < 2e-5, "{text}");
    assert!((field("mean_jump_k") + 0.5).abs() < 2e-6, "{text}");
    assert!(field("relative_closure") < 1e-6, "{text}");

    // f85xj.6.8: `report` exports the published field as VTU, byte-for-byte
    // the artifact the conduction receipt names; the independent checker
    // reads it and its extrema and cell count are the receipt's own.
    let exported = run(args(&[
        "--json",
        "report",
        run_id,
        dir.join("pair.db").to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        exported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        exported.stderr
    );
    assert!(
        exported.stdout.contains("\"field_vtu\":"),
        "{}",
        exported.stdout
    );
    let vtu = std::fs::read(dir.join(format!("{run_id}.field.vtu"))).expect("field exported");
    let named = text
        .split("\"field_artifact\":\"")
        .nth(1)
        .map(|rest| &rest[..64])
        .expect("receipt names the field");
    assert_eq!(fs_blake3::hash_bytes(&vtu).to_hex(), named);
    let checked =
        fs_viz::vtu::VtuChecker::check(std::str::from_utf8(&vtu).unwrap()).expect("VTU checks");
    assert_eq!(checked.num_cells as f64, field("elements"));
    let (_, [t_lo, t_hi]) = checked
        .array_extrema
        .iter()
        .find(|(name, _)| name == "temperature")
        .expect("temperature array")
        .clone();
    let temperature = text
        .split("\"temperature\":{\"unit\":\"K\",")
        .nth(1)
        .expect("receipt temperature object");
    assert_eq!(t_lo, number_after(temperature, "\"min\":"));
    assert_eq!(t_hi, number_after(temperature, "\"max\":"));
}

#[test]
fn g1_compare_answers_the_heatsink_fan_speed_decision_with_pressure_drop() {
    // q61wp.83 first slice: the heatsink example's declared decision is
    // "compare fan operating points". It now requests `pressure-drop` beside
    // `temperature-max`, and `compare` diffs both. The exact falsifier is the fan
    // affinity law: against the quadratic orifice/leakage network, the
    // operating pressure scales with speed squared, so 0.7 -> 0.9 must raise the
    // pressure drop by (0.9/0.7)^2 while the faster air cools the part.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("fan-speed-decision");
    let ledger = dir.join("fans.db");
    let fast = dir.join("heatsink-fan-0p9.fsim");
    let declared = std::fs::read_to_string(&fsim).unwrap();
    assert_eq!(declared.matches(":speed-ratio 0.7").count(), 1);
    std::fs::write(
        &fast,
        declared.replace(":speed-ratio 0.7", ":speed-ratio 0.9"),
    )
    .unwrap();
    let mut runs = Vec::new();
    for project in [&fsim, &fast] {
        let imported = run(args(&[
            "--json",
            "import",
            project.to_string_lossy().as_ref(),
            stl.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]));
        assert_eq!(
            imported.exit_code,
            exit::SUCCESS,
            "stderr: {}",
            imported.stderr
        );
        let solved = run(args(&[
            "--json",
            "solve",
            project.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--materials",
            pack.to_string_lossy().as_ref(),
        ]));
        assert_eq!(solved.exit_code, exit::SUCCESS, "stderr: {}", solved.stderr);
        runs.push(
            solved
                .stdout
                .split("\"run\":\"")
                .nth(1)
                .and_then(|rest| rest.get(..64))
                .unwrap()
                .to_string(),
        );
    }
    let compared = run(args(&[
        "--json",
        "compare",
        &runs[0],
        &runs[1],
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        compared.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        compared.stderr
    );
    let diff = |name: &str| -> (f64, f64) {
        let row = compared
            .stdout
            .split(&format!("\"name\":\"{name}\""))
            .nth(1)
            .unwrap_or_else(|| panic!("no {name} diff in {}", compared.stdout));
        (
            number_after(row, "\"nominal_left\":"),
            number_after(row, "\"nominal_right\":"),
        )
    };
    let (dp_slow, dp_fast) = diff("pressure-drop");
    let ratio = 0.9_f64 / 0.7;
    let affinity = ratio * ratio;
    assert!(
        (dp_fast / dp_slow / affinity - 1.0).abs() < 1e-9,
        "pressure drop {dp_slow} -> {dp_fast} Pa is not the affinity-law ratio {affinity}"
    );
    let (t_slow, t_fast) = diff("temperature-max");
    assert!(
        t_fast < t_slow,
        "faster air must cool the part: {t_slow} -> {t_fast} K"
    );
    // Fan input power dp * Q / eta scales with speed cubed at fixed efficiency.
    let (p_slow, p_fast) = diff("fan-power");
    assert!(
        (p_fast / p_slow / (affinity * ratio) - 1.0).abs() < 1e-9,
        "fan power {p_slow} -> {p_fast} W is not the cube-law ratio"
    );

    // Fan power without a cited efficiency refuses by name; it never
    // assumes 100 % or any default.
    let start = declared
        .find("(fan-efficiency ")
        .expect("example cites a fan efficiency");
    let end = start + declared[start..].find(')').unwrap() + 1;
    let uncited = dir.join("heatsink-fan-uncited.fsim");
    std::fs::write(
        &uncited,
        format!("{}{}", declared[..start].trim_end(), &declared[end..]),
    )
    .unwrap();
    let imported = run(args(&[
        "--json",
        "import",
        uncited.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );
    let refused = run(args(&[
        "--json",
        "solve",
        uncited.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert!(
        refused
            .stderr
            .contains("cli-solve-qoi-fan-power-no-efficiency"),
        "stderr: {}",
        refused.stderr
    );
}

#[test]
fn g1_surface_mean_over_the_whole_skin_equals_the_energy_balance_exactly() {
    // q61wp.83 slice 3: the surface family over a declared surface. The
    // plate-hole example declares a `skin` surface covering every exterior
    // face (hole walls included) and requests its area mean. All 2 W leave
    // by uniform h = 10 W/m^2/K convection to 293.15 K, so the area mean is
    // EXACTLY 293.15 + P/(h A) with A = 5504 mm^2 from the generator. The P1
    // face-integral mean and the Robin flux integral see the same field, so
    // only the solver tolerance separates them.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/plate-hole/plate-hole.fsim");
    let stl = root.join("examples/plate-hole/plate-hole.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("plate-hole-skin");
    let ledger = dir.join("plate.db");
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );
    let solved = run(args(&[
        "--json",
        "solve",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(solved.exit_code, exit::SUCCESS, "stderr: {}", solved.stderr);
    let run_id = solved
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.get(..64))
        .unwrap();
    let ledger = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap();
    let receipts = stage_receipt_hashes(&ledger, run_id);
    let qoi = receipt_text(&ledger, &receipts[5]);
    let row = qoi
        .split("\"name\":\"surface-mean-temperature\"")
        .nth(1)
        .unwrap_or_else(|| panic!("no surface-mean row in {qoi}"));
    assert!(row.contains("\"region\":\"skin\""), "{row}");
    let mean = number_after(row, "\"value\":");
    let exact = 293.15 + 2.0 / (10.0 * 5.504e-3);
    // TOLERANCE 5e-6 K: MEASURED 2026-09-30 diff 1.05e-6 K (solver
    // tolerance-rel 1e-6 on a 36 K rise), ~5x headroom.
    eprintln!(
        "skin mean {mean} K vs exact {exact} K (diff {:e})",
        mean - exact
    );
    assert!(
        (mean - exact).abs() < 5e-6,
        "skin mean {mean} K is not the energy balance {exact} K"
    );
}

#[test]
fn g1_a_passive_heatsink_converges_on_the_natural_convection_card() {
    // fsim v10: the same heatsink with no fan, cooled by buoyant air through
    // the Churchill-Chu vertical-plate card. The coefficient depends on the
    // solved wall-to-ambient difference, so the stage iterates to a fixed
    // point. Independent checks: all declared power leaves through the law;
    // the receipt's Nu is the Churchill-Chu formula, written out here
    // independently, at the receipt's own Ra; h = Nu k / L; and doubling the
    // power raises the difference by LESS than 2x (h grows with it).
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-natural.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("heatsink-natural");
    let solve = |project: &std::path::Path, tag: &str| -> String {
        let ledger = dir.join(format!("{tag}.db"));
        let imported = run(args(&[
            "--json",
            "import",
            project.to_string_lossy().as_ref(),
            stl.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ]));
        assert_eq!(
            imported.exit_code,
            exit::SUCCESS,
            "stderr: {}",
            imported.stderr
        );
        let solved = run(args(&[
            "--json",
            "solve",
            project.to_string_lossy().as_ref(),
            ledger.to_string_lossy().as_ref(),
            "--materials",
            pack.to_string_lossy().as_ref(),
        ]));
        assert_eq!(solved.exit_code, exit::SUCCESS, "stderr: {}", solved.stderr);
        let run_id = solved
            .stdout
            .split("\"run\":\"")
            .nth(1)
            .and_then(|rest| rest.get(..64))
            .unwrap()
            .to_string();
        let ledger = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap();
        let receipts = stage_receipt_hashes(&ledger, &run_id);
        receipt_text(&ledger, &receipts[4])
    };
    let text = solve(&fsim, "three");
    let law = text
        .split("\"natural_convection\":")
        .nth(1)
        .unwrap_or_else(|| panic!("no natural_convection block in {text}"));
    let field = |key: &str| number_after(law, &format!("\"{key}\":"));
    let (htc, delta_t, rayleigh, nusselt) = (
        field("htc_w_m2_k"),
        field("delta_t_k"),
        field("rayleigh"),
        field("nusselt"),
    );
    eprintln!(
        "natural: h {htc} W/m2K, dT {delta_t} K, Ra {rayleigh:e}, Nu {nusselt}, iterations {}",
        field("iterations")
    );
    assert!(
        (field("heat_rate_w") - 3.0).abs() < 1e-6,
        "all 3 W must leave by natural convection: {law}"
    );
    // Churchill & Chu (1975), full-range vertical plate, Pr = 0.707.
    let pr = 0.707_f64;
    let shape = (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let root_nu = 0.825 + 0.387 * rayleigh.powf(1.0 / 6.0) / shape;
    let churchill_chu = root_nu * root_nu;
    assert!(
        (nusselt / churchill_chu - 1.0).abs() < 1e-9,
        "Nu {nusselt} vs Churchill-Chu {churchill_chu}"
    );
    assert!(
        (htc / (nusselt * 26.3e-3 / 0.06) - 1.0).abs() < 1e-12,
        "h = Nu k / L"
    );
    assert!(delta_t > 0.0 && field("iterations") < 80.0);

    let declared = std::fs::read_to_string(&fsim).unwrap();
    assert_eq!(declared.matches(":watts 3.0kg").count(), 1);
    let doubled = dir.join("heatsink-natural-6w.fsim");
    std::fs::write(&doubled, declared.replace(":watts 3.0kg", ":watts 6.0kg")).unwrap();
    let text = solve(&doubled, "six");
    let law = text.split("\"natural_convection\":").nth(1).unwrap();
    let delta_6 = number_after(law, "\"delta_t_k\":");
    let ratio = delta_6 / delta_t;
    eprintln!("natural: dT 3 W {delta_t} K, 6 W {delta_6} K, ratio {ratio}");
    assert!(
        ratio > 1.0 && ratio < 2.0,
        "doubling power must raise dT by less than 2x: {ratio}"
    );
}

#[test]
fn g1_a_passive_heatsink_also_radiates_and_runs_cooler() {
    // A fanless heatsink sheds a comparable share by radiation. The gray
    // reference surface on the metal augments the natural-convection row: the
    // fixed point must still converge on the CONVECTIVE part, the convective
    // and radiative watts must account for all 3 W, and the part must run
    // cooler than with natural convection alone.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let natural =
        std::fs::read_to_string(root.join("examples/heatsink-fan/heatsink-natural.fsim")).unwrap();
    let law_end = ":correlation \"convection.churchill-chu-vertical-plate\"))";
    assert_eq!(natural.matches(law_end).count(), 1);
    let radiating = natural.replace(
        law_end,
        ":correlation \"convection.churchill-chu-vertical-plate\")) :radiation (radiation :surfaces (surfaces (surface :name \"gray-metal\" :target \"metal\" :card \"63485429663ba24d53d67a7d3b03ab0611f17f1e3e445e6a1ef4f636093a6e4f\" :query-temperature 300.0K :reservoir-temperature 293.15K)) :max-iterations 128 :temperature-tolerance 1e-8K :heat-tolerance 1e-7kg·m^2·s^-3 :relaxation 0.5)",
    );
    let dir = scratch("heatsink-natural-radiating");
    let fsim = dir.join("heatsink-natural-radiating.fsim");
    std::fs::write(&fsim, &radiating).unwrap();
    let ledger = dir.join("rad.db");
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        root.join("examples/heatsink-fan/heatsink.stl")
            .to_string_lossy()
            .as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );
    let solved = run(args(&[
        "--json",
        "solve",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        root.join("data/reference-project/aa6061.fsmcdpk")
            .to_string_lossy()
            .as_ref(),
        "--materials",
        root.join("data/reference-project/gray-surface.fsmcdpk")
            .to_string_lossy()
            .as_ref(),
    ]));
    assert_eq!(
        solved.exit_code,
        exit::SUCCESS,
        "stdout {} stderr {}",
        solved.stdout,
        solved.stderr
    );
    let run_id = solved
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.get(..64))
        .unwrap();
    let ledger = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap();
    let receipts = stage_receipt_hashes(&ledger, run_id);
    let text = receipt_text(&ledger, &receipts[4]);
    let law = text
        .split("\"natural_convection\":")
        .nth(1)
        .expect("natural block");
    let convective = number_after(law, "\"heat_rate_w\":");
    let radiative = number_after(&text, "\"radiative_out_w\":");
    let maximum = number_after(&text, "\"max\":");
    eprintln!(
        "natural+radiation: convective {convective} W, radiative {radiative} W, T_max {maximum} K"
    );
    assert!(radiative > 0.0 && convective > 0.0, "{text}");
    assert!(
        (convective + radiative - 3.0).abs() < 1e-5,
        "the two exits carry all 3 W: {convective} + {radiative}"
    );
    assert!(
        maximum < 316.73,
        "radiation must cool the passive part below 316.74 K, got {maximum}"
    );
}

#[test]
fn g0_package_missing_ledger_fails_closed() {
    let output = run(args(&[
        "package",
        "0000000000000000000000000000000000000000000000000000000000000000",
        "/nonexistent/ledger.db",
        "--json",
    ]));
    assert_eq!(output.exit_code, exit::INPUT, "stderr: {}", output.stderr);
    assert!(output.stderr.contains("cli-export-ledger-missing"));
    assert!(!output.stdout.contains("\"verdict\":\"pass\""));
    assert!(
        !std::path::Path::new("/nonexistent/ledger.db").exists(),
        "an export must never create a ledger"
    );
}

#[test]
fn g0_empty_ledger_cannot_mint_a_self_consistent_package() {
    let dir = scratch("package");
    let ledger = dir.join("test_ledger.db");
    let _ = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap();

    let run_id = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
    let output = run(args(&[
        "package",
        run_id,
        ledger.to_string_lossy().as_ref(),
        "--json",
    ]));
    assert_eq!(output.exit_code, exit::REFUSED, "stderr: {}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"refused\""));
    assert!(output.stderr.contains("cli-solve-unknown-run"));
    assert!(!output.stdout.contains("\"merkle_root\""));
    assert!(!output.stdout.contains("\"verdict\":\"pass\""));
    assert!(!output.stdout.contains("junction_maximum"));
    assert!(!dir.join(format!("{run_id}.fspkg")).exists());
}

#[test]
fn g0_report_missing_ledger_fails_closed() {
    let output = run(args(&[
        "report",
        "0000000000000000000000000000000000000000000000000000000000000000",
        "/nonexistent/ledger.db",
        "--json",
    ]));
    assert_eq!(output.exit_code, exit::INPUT, "stderr: {}", output.stderr);
    assert!(output.stderr.contains("cli-export-ledger-missing"));
    assert!(!output.stdout.contains("junction_maximum"));
}

#[test]
fn g0_empty_ledger_cannot_mint_a_verified_engineering_report() {
    let dir = scratch("report");
    let ledger = dir.join("report_test_ledger.db");
    let _ = fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap();

    let run_id = "feedface000000000000000000000000feedface000000000000000000000000";
    let output = run(args(&[
        "report",
        run_id,
        ledger.to_string_lossy().as_ref(),
        "--json",
    ]));
    assert_eq!(output.exit_code, exit::REFUSED, "stderr: {}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"refused\""));
    assert!(output.stderr.contains("cli-solve-unknown-run"));
    assert!(!output.stdout.contains("\"content_hash\""));
    assert!(!output.stdout.contains("junction_maximum"));
    assert!(!output.stdout.contains("Verified"));

    let html_path = dir.join(format!("{run_id}.report.html"));
    let json_path = dir.join(format!("{run_id}.report.json"));

    assert!(!html_path.exists(), "a refused report must not write HTML");
    assert!(
        !json_path.exists(),
        "a refused report must not write a JSON twin"
    );
}

#[test]
fn g1_run_completes_seven_stages_and_exports_report_and_package_for_the_reference_project() {
    // The tracked reference project declares conduction and a temperature
    // maximum, so the real binary must now carry it through all seven solve
    // stages, seal the report stage in the ledger, and export the retained
    // report, JSON twin, and evidence package. Every displayed value traces to
    // a retained receipt; nothing here is allowed to be a literal.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("data/reference-project/cooling-reference.fsim");
    let stl = root.join("data/reference-project/plate.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("run-complete");
    let ledger = dir.join("complete.db");

    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );

    let output = run(args(&[
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    assert!(output.stdout.contains("\"command\":\"run\""));
    assert!(output.stdout.contains("\"status\":\"completed\""));
    assert!(output.stdout.contains("\"stages_completed\":7"));
    assert!(output.stdout.contains("\"checker\":\"pass\""));
    assert!(
        output
            .stderr
            .contains("\"stage\":\"report\",\"ordinal\":6,\"status\":\"ok\""),
        "the report stage reports progress like every other stage: {}",
        output.stderr
    );
    let run_id = output
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .expect("run id in the result")
        .to_string();
    assert_eq!(run_id.len(), 64);

    let html = std::fs::read_to_string(dir.join(format!("{run_id}.report.html")))
        .expect("the retained HTML report was exported");
    let twin = std::fs::read_to_string(dir.join(format!("{run_id}.report.json")))
        .expect("the retained JSON twin was exported");
    let package = std::fs::read_to_string(dir.join(format!("{run_id}.fspkg")))
        .expect("the retained package was exported");
    assert!(html.contains(&run_id));
    assert!(html.contains("temperature-max"));
    assert!(
        html.contains("NO-DATA"),
        "unmeasured terms print as NO-DATA, never as numbers"
    );
    assert!(html.contains("Estimated"));
    assert!(!html.contains("342.15"));
    assert!(twin.contains("\"schema\": \"frankensim.report.engineering.v1\""));
    assert!(twin.contains("\"state\": \"no-data\""));
    assert!(twin.contains("\"stage\": \"qoi\""));
    assert!(!twin.contains("NaN"), "the JSON twin never emits NaN");
    let parsed = fs_package::EvidencePackage::from_json(&package).expect("format-9 package");
    assert!(fs_checker::check(&parsed).passed());

    // Exports are idempotent: identical bytes already on disk are accepted.
    let again = run(args(&[
        "--json",
        "report",
        &run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(again.exit_code, exit::SUCCESS, "stderr: {}", again.stderr);
    assert!(again.stdout.contains("\"stages_completed\":7"));
    assert!(
        again
            .stdout
            .contains("\"verification\":\"sealed-evidence\""),
        "exports prove the run by sealed evidence, never by replaying physics: {}",
        again.stdout
    );
    let packaged = run(args(&[
        "--json",
        "package",
        &run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        packaged.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        packaged.stderr
    );
    assert!(packaged.stdout.contains("\"checker\":\"pass\""));
    assert!(packaged.stdout.contains("\"merkle_root\":\""));

    // A differing file at the export path is a conflict, never overwritten.
    std::fs::write(dir.join(format!("{run_id}.fspkg")), b"tampered").expect("tamper");
    let conflict = run(args(&[
        "--json",
        "package",
        &run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(conflict.exit_code, exit::REFUSED);
    assert!(conflict.stderr.contains("cli-export-output-conflict"));
    assert_eq!(
        std::fs::read(dir.join(format!("{run_id}.fspkg"))).expect("still there"),
        b"tampered"
    );
}

#[test]
fn g3_report_json_conflict_does_not_publish_a_partial_html_twin() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("data/reference-project/cooling-reference.fsim");
    let stl = root.join("data/reference-project/plate.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("report-twin-conflict");
    let ledger = dir.join("report_twin_conflict.db");

    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(imported.exit_code, exit::SUCCESS, "{}", imported.stderr);

    // `solve` seals all seven stages without running the workflow exports, so
    // both twin destinations begin absent.
    let solved = run(args(&[
        "--json",
        "solve",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(solved.exit_code, exit::SUCCESS, "{}", solved.stderr);
    let run_id = solved
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .expect("completed solve reports its run identity");
    let html_path = dir.join(format!("{run_id}.report.html"));
    let json_path = dir.join(format!("{run_id}.report.json"));
    assert!(!html_path.exists());
    assert!(!json_path.exists());

    std::fs::write(&json_path, b"conflicting JSON twin").expect("hostile twin writes");
    let refused = run(args(&[
        "--json",
        "report",
        run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(refused.exit_code, exit::REFUSED, "{}", refused.stderr);
    assert!(refused.stderr.contains("cli-export-output-conflict"));
    assert!(
        !html_path.exists(),
        "a known JSON conflict must refuse before the HTML twin is published"
    );
    assert_eq!(
        std::fs::read(json_path).expect("hostile twin remains"),
        b"conflicting JSON twin"
    );
}

#[test]
fn g0_run_stops_at_the_conduction_gap_when_the_project_declares_no_conduction() {
    // Strip the conduction declaration from the heatsink example in a scratch
    // copy: `run` must then refuse at the conduction stage by name (exit 4,
    // `cli-solve-conduction-undeclared` — a project defect, not a stage gap),
    // name the stage, and write no report or package. This pins the negative
    // space of the example above: conduction executes only for a declared
    // solid problem, never by inventing one.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let source = std::fs::read_to_string(root.join("examples/heatsink-fan/heatsink-fan.fsim"))
        .expect("example is readable");
    let start = source
        .find("(conduction ")
        .expect("example declares conduction");
    let mut depth = 0usize;
    let mut end = None;
    for (offset, ch) in source[start..].char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    end = Some(start + offset + 1);
                    break;
                }
            }
            _ => {}
        }
    }
    let end = end.expect("balanced conduction form");
    let stripped = format!("{}{}", &source[..start], &source[end..])
        .replace(" )", ")")
        .replace("  (", " (");
    let dir = scratch("run-no-conduction");
    let fsim = dir.join("heatsink-no-conduction.fsim");
    std::fs::write(&fsim, stripped.trim_end()).expect("scratch project");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let ledger = dir.join("no-conduction.db");

    let validated = run(args(&[
        "--json",
        "validate",
        fsim.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        validated.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        validated.stderr
    );
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );

    let output = run(args(&[
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
        "--json",
    ]));
    assert_eq!(
        output.exit_code,
        exit::REFUSED,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    assert!(
        output.stderr.contains("cli-solve-conduction-undeclared"),
        "stderr: {}",
        output.stderr
    );
    assert!(
        output.stdout.contains("\"stage\":\"conduction\""),
        "stdout: {}",
        output.stdout
    );
    assert!(!output.stdout.contains("\"status\":\"completed\""));
    assert!(!output.stdout.contains("\"report_html\""));
    let exported: Vec<_> = std::fs::read_dir(&dir)
        .expect("scratch dir")
        .filter_map(Result::ok)
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|name| name.ends_with(".report.html") || name.ends_with(".fspkg"))
        .collect();
    assert!(
        exported.is_empty(),
        "a gapped run must export nothing: {exported:?}"
    );
}

/// Hex content hashes of every `solve-stage-receipt` the run retained, in op
/// order (import-verify, assign, material-resolve, flow-network, conduction,
/// qoi, report).
fn stage_receipt_hashes(ledger: &fs_ledger::Ledger, run_hex: &str) -> Vec<String> {
    let run = fs_cli::SolveRunId::parse_hex(run_hex).expect("hex run id");
    let mut ids = ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)
        .expect("ops");
    ids.sort_unstable();
    let mut receipts = Vec::new();
    for id in ids {
        let Some(row) = ledger.op(id).expect("op row") else {
            continue;
        };
        if row.session.as_deref() != Some(run.as_bytes().as_slice())
            || row.outcome.as_deref() != Some("ok")
        {
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

fn receipt_text(ledger: &fs_ledger::Ledger, hex: &str) -> String {
    let hash = fs_ledger::ContentHash::from_hex(hex).expect("content hash");
    let bytes = ledger
        .get_artifact(&hash)
        .expect("read")
        .expect("retained receipt");
    String::from_utf8(bytes).expect("receipt is UTF-8")
}

/// The number that follows `marker` in a receipt.
fn number_after(text: &str, marker: &str) -> f64 {
    text.split(marker)
        .nth(1)
        .and_then(|rest| rest.split([',', '}']).next())
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or_else(|| panic!("numeric field after `{marker}` in {text}"))
}

/// Journey A physical anchor (bead frankensim-rc-root-q61wp.14 item 4): a
/// Level-A hand calculation from the same project and card, with no
/// independent loose tolerance.
///
/// The reference plate dissipates `Q` uniformly and loses it through one
/// uniform convection boundary `h` to `T_ref`, with no other heat path. At
/// steady state the surface energy balance is an identity, not a model:
/// `h · ∮ (T − T_ref) dA = Q`, so the area-weighted mean surface excess is
/// exactly `Q / (h · A)` and the solver's minimum and maximum must bracket
/// `T_ref + Q / (h · A)`. That bracket has no tolerance. The width of the
/// bracket is the solid's internal spread, which conduction bounds by the
/// Biot number `h · L / k`: with aluminium (`k` read from the card) and
/// `h = 10 W/m²K` the spread must be a small fraction of the lumped excess.
/// The receipt's own energy block must close the same balance. Every input
/// is read from the project, the card, and the STL; every output from the
/// retained conduction and QoI receipts.
#[test]
fn ja_005_level_a_energy_balance_brackets_the_retained_maximum() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("data/reference-project/cooling-reference.fsim");
    let stl = root.join("data/reference-project/plate.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch("ja-005");
    let ledger_path = dir.join("anchor.db");

    // Inputs from the project itself.
    let source = std::fs::read_to_string(&fsim).expect("reference project reads");
    let decoded = fs_project::parse_sexpr_migrating(&source)
        .expect("historical reference project migrates")
        .decoded;
    let power = decoded.spec.power.as_ref().expect("power declared");
    assert_eq!(power.len(), 1, "one dissipating region");
    let q_w = power[0].watts.value * power[0].duty;
    let conduction = decoded
        .spec
        .cooling
        .as_ref()
        .expect("cooling declared")
        .conduction
        .as_ref()
        .expect("conduction declared");
    assert_eq!(conduction.boundaries.len(), 1, "one thermal boundary");
    assert!(!conduction.adiabatic_remainder);
    let (h, t_ref) = match &conduction.boundaries[0].condition {
        fs_project::spec::ThermalBoundaryCondition::Convection {
            coefficient,
            reference_temperature,
        } => (coefficient.value, reference_temperature.value),
        other => {
            panic!("the reference plate declares a coefficient convection law, found {other:?}")
        }
    };
    assert!(q_w > 0.0 && h > 0.0 && t_ref > 0.0);

    // Conductivity from the card the project binds.
    let card = fs_matdb::NormalizedMaterialCardPack::from_bytes(
        &std::fs::read(&pack).expect("card pack reads"),
    )
    .expect("card pack parses");
    let claims = card.card().claims_for("thermal-conductivity");
    assert_eq!(claims.len(), 1, "one conductivity claim");
    let k = match &claims[0].1.value {
        fs_matdb::PropertyValue::Scalar { value, .. } => *value,
        other => panic!("scalar conductivity expected, found {other:?}"),
    };
    assert!(k > 0.0);

    // Wetted area and characteristic length from the STL.
    let soup = fs_io::quarantine::import_mesh(&std::fs::read(&stl).expect("stl reads"), "stl")
        .expect("plate imports")
        .into_inner();
    let mut area = 0.0_f64;
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for point in &soup.positions {
        for (axis, value) in [point.x, point.y, point.z].into_iter().enumerate() {
            lo[axis] = lo[axis].min(value);
            hi[axis] = hi[axis].max(value);
        }
    }
    for t in 0..soup.triangles.len() {
        let [a, b, c] = soup.tri(t);
        let u = [b.x - a.x, b.y - a.y, b.z - a.z];
        let v = [c.x - a.x, c.y - a.y, c.z - a.z];
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        area += 0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
    }
    let length = (0..3)
        .map(|axis| hi[axis] - lo[axis])
        .fold(0.0_f64, f64::max);
    assert!(area > 0.0 && length > 0.0);
    let biot = h * length / k;
    let lumped_excess_k = q_w / (h * area);

    // The solver's retained answer.
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(imported.exit_code, exit::SUCCESS, "{}", imported.stderr);
    let output = run(args(&[
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    let run_id = output
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .expect("run id")
        .to_string();
    let ledger =
        fs_ledger::Ledger::open(ledger_path.to_str().expect("utf-8 path")).expect("ledger");
    let receipts = stage_receipt_hashes(&ledger, &run_id);
    assert_eq!(receipts.len(), 7, "seven stage receipts");
    let conduction_receipt = receipt_text(&ledger, &receipts[4]);
    assert!(conduction_receipt.contains("\"stage\":\"conduction\""));
    let qoi_receipt = receipt_text(&ledger, &receipts[5]);
    assert!(qoi_receipt.contains("\"stage\":\"qoi\""));
    let t_min = number_after(
        &conduction_receipt,
        "\"temperature\":{\"unit\":\"K\",\"min\":",
    );
    let t_max = number_after(&conduction_receipt, "\"max\":");
    let source_w = number_after(&conduction_receipt, "\"source_w\":");
    let robin_out_w = number_after(&conduction_receipt, "\"robin_out_w\":");
    let closure_w = number_after(&conduction_receipt, "\"closure_w\":");
    let qoi_value = number_after(&qoi_receipt, "\"value\":");

    // (1) The receipt's energy block is the same balance the hand calc uses.
    assert!(
        (source_w - q_w).abs() <= 1e-9 * q_w,
        "retained source {source_w} W is the declared duty {q_w} W"
    );
    assert!(
        closure_w.abs() <= 1e-6 * q_w,
        "the retained balance closes: closure {closure_w} W of {q_w} W"
    );
    assert!(
        (robin_out_w - q_w).abs() <= 1e-6 * q_w,
        "all heat leaves through the convection boundary: {robin_out_w} W of {q_w} W"
    );

    // (2) The exact bracket: min excess <= Q/(hA) <= max excess.
    let min_excess = t_min - t_ref;
    let max_excess = t_max - t_ref;
    assert!(
        min_excess <= lumped_excess_k && lumped_excess_k <= max_excess,
        "the solver's surface temperatures must bracket the lumped excess: \
         min {min_excess} K <= Q/(hA) = {lumped_excess_k} K <= max {max_excess} K"
    );
    assert_eq!(qoi_value, t_max, "the QoI is the retained maximum");

    // (3) The bracket width is the internal conduction spread, bounded by
    // the Biot number: for Bi << 1 the body is near-isothermal, so the
    // spread is a small fraction of the lumped excess. The factor is the
    // order-one geometry constant of a 1-D slab with distributed source
    // (spread <= Q L / (2 k A) = Bi/2 · Q/(hA)); assert it with 2x headroom.
    let spread = t_max - t_min;
    assert!(
        biot < 0.5,
        "the anchor is a lumped calculation; Bi = {biot} must be small"
    );
    assert!(
        spread <= biot * lumped_excess_k,
        "internal spread {spread} K exceeds the Biot bound {} K (Bi = {biot})",
        biot * lumped_excess_k
    );
    println!(
        "{{\"falsifier\":\"level-a-energy-balance-anchor\",\"q_w\":{q_w},\"h_w_m2k\":{h},\"t_ref_k\":{t_ref},\"k_w_mk\":{k},\"area_m2\":{area},\"length_m\":{length},\"biot\":{biot},\"lumped_excess_k\":{lumped_excess_k},\"t_min_k\":{t_min},\"t_max_k\":{t_max},\"spread_k\":{spread},\"closure_w\":{closure_w}}}"
    );
}
