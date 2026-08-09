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
            .contains("import <project> <source> <ledger.db>")
    );

    let duplicate = run(args(&["validate", "x.fsim", "--json", "--json"]));
    assert_eq!(duplicate.exit_code, exit::USAGE);
    assert!(duplicate.stderr.contains("cli-duplicate-flag"));

    let extra = run(args(&["report", "run-1", "extra"]));
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
fn g0_unintegrated_product_stages_fail_closed_with_their_owner() {
    for (command, dependency) in [
        (&["report", "run-1"][..], "f85xj.6.9"),
        (&["package", "run-1"][..], "f85xj.6.10"),
    ] {
        let mut invocation = args(command);
        invocation.push("--json".to_string());
        let output = run(invocation);
        assert_eq!(output.exit_code, exit::UNAVAILABLE, "{command:?}");
        assert!(output.stdout.contains("\"status\":\"unavailable\""));
        assert!(output.stdout.contains(dependency), "{command:?}");
        assert!(output.stderr.contains("cli-stage-unavailable"));
        assert!(output.stderr.contains("placeholder artifact"));
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
    // Pinned, not endorsed: this is a resource cap reported in the USAGE exit
    // class even though the invocation matches the documented grammar and
    // `exit::INPUT` is the class documented for "admitted by the CLI resource
    // cap". Tracked as its own bead; this assertion exists so the class cannot
    // drift silently while that is decided.
    assert_eq!(past_cap.exit_code, exit::USAGE);
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
