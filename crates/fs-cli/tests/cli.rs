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
            conduction: None,
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
    // Frozen canonical hashes: any fixture-byte drift fails right here
    // instead of rotting quietly. Regenerate with the real verb and commit
    // fixture + hash in the same commit.
    assert!(
        output.stdout.contains(
            "\"project_hash\":\"1f135da400dec2bed1bba833b033026ce37bc93c113c79c25f6c6fb4e730780b\""
        ),
        "heated-plate.fsim drifted from its frozen canonical hash"
    );

    let reference = root.join("data/reference-project/cooling-reference.fsim");
    let ref_out = run(args(&["--json", "validate", &reference.to_string_lossy()]));
    assert_eq!(ref_out.exit_code, exit::SUCCESS);
    assert!(
        ref_out.stdout.contains(
            "\"project_hash\":\"4e2d71ab877ec805b7aa617d5d8bd2ca70f6ea38ca5afd8f0742b23ae60d7135\""
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
        again.stdout.contains("\"verification\":\"sealed-evidence\""),
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
    let stripped = format!("{}{}", &source[..start], &source[end..]).replace(" )", ")");
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
