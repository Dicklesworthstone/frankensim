//! G0 command-contract tests for the `frankensim` CLI membrane.

use fs_cli::{exit, run, validate_source};
use fs_project::{
    Budgets, Cooling, EntityDecl, Envelope, GeometryArtifact, GeometryAssignment, MeshSelector,
    Metadata, OutputRequest, PowerDissipation, ProjectSpec, Seeds, SolverSettings, ThermalLimit,
    UnitsDoctrine, Versions, print_sexpr,
};
use fs_qty::QtyAny;

fn args(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
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
        }),
        versions: Some(Versions {
            schema: fs_project::FSIM_VERSION,
            constellation: "00".repeat(32),
            workspace: "11".repeat(20),
        }),
        seeds: Some(Seeds { master: 7 }),
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
        power: Some(vec![PowerDissipation {
            region: "hot".to_string(),
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
            class: "surface".to_string(),
            region: "hot".to_string(),
            limit: kelvin(353.15),
            margin: kelvin(5.0),
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
        (&["solve", "project.fsim"][..], "f85xj.6.5"),
        (&["solve", "--resume", "run-1"][..], "f85xj.6.5"),
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
