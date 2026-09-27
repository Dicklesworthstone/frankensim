//! Fan uncertainty reaches the existing flow-network and conjugate solid/air
//! solve. Replay keeps the original curve, bank domain and absolute ratio.

use super::*;

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    let cooling = fixture.project.cooling.as_mut().unwrap();
    // A nonunit baseline distinguishes setting the source-curve ratio from
    // accidentally multiplying the base speed for each sample.
    cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio = 0.7;
    cooling.conduction.as_mut().unwrap().boundaries[0].condition =
        fs_project::ThermalBoundaryCondition::AirflowConvection {
            branch: "air".into(),
            order: 0,
            inlet_temperature: fs_qty::QtyAny::new(294.0, fs_project::spec::dims::TEMPERATURE),
            hydraulic_diameter: fs_qty::QtyAny::new(0.02, fs_project::spec::dims::LENGTH),
            flow_area: fs_qty::QtyAny::new(0.004, fs_project::spec::dims::AREA),
            channel_length: fs_qty::QtyAny::new(0.3, fs_project::spec::dims::LENGTH),
            correlation: "convection.gnielinski".into(),
        };
    std::fs::write(
        fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap(),
    ).unwrap();
    let study = STUDY.replace(
        "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
        "(uniform :name \"speed\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.8 :high 1.2)",
    ).replace(
        ":target convection-temperature :entity \"air\" :low 294K :high 300K",
        ":target air-inlet-temperature :entity \"air\" :low 294K :high 294K",
    );
    std::fs::write(fixture.sources.join("study.fsim"), study).unwrap();
    fixture
}

#[test]
fn g1_native_fan_speed_samples_match_direct_coupled_solves_and_change_temperature() {
    let fixture = fixture();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    assert_eq!(result.str_field("status"), Some("completed"));
    let mut responses = Vec::new();
    for (ordinal, row) in rows(&report).iter().enumerate() {
        let parameters = row.get("parameters").unwrap().as_array().unwrap();
        let speed = parameters[0].as_f64().unwrap();
        assert!((0.8..=1.2).contains(&speed));
        assert_eq!(parameters[1].as_f64(), Some(294.0));
        // Bind the independent native solve by hand, bypassing the uncertainty
        // parser's target application and sample-project producer.
        let mut project = fixture.project.clone();
        project.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = speed;
        let (project_hash, run, temperature) = fixture.solve_project(ordinal, &project);
        assert_eq!(row.str_field("project_hash"), Some(project_hash.as_str()));
        assert_eq!(row.str_field("run"), Some(run.as_str()));
        assert_eq!(row.f64_field("value_k").unwrap().to_bits(), temperature.to_bits());
        responses.push((speed, temperature));
    }
    responses.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert!(responses.last().unwrap().0 - responses[0].0 > 0.01);
    for pair in responses.windows(2) {
        assert!(pair[0].1 > pair[1].1,
            "higher sampled speed must cool this fixed-power conjugate fixture: {pair:?}");
    }
    assert!(responses[0].1 - responses.last().unwrap().1 > 1e-8);
}

#[test]
fn g4_g5_native_fan_speed_resume_retains_original_bank_and_replays_exactly() {
    let whole = fixture();
    let result = whole.study(None, fs_cli::exit::SUCCESS);
    let (expected, _) = whole.report(&result);
    let split = fixture();
    let partial = split.study(Some("2"), fs_cli::exit::BUDGET);
    let (_, before) = split.report(&partial);
    assert_eq!(rows(&before).len(), 2);
    std::fs::rename(&split.sources, split.dir.join("relocated-fan-sources")).unwrap();
    let result = split.resume(partial.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual, report) = split.report(&result);
    assert_eq!(actual, expected, "resumption uses the retained curve, bank and sampling ordinals");
    assert_eq!(&rows(&report)[..2], rows(&before));
}

#[test]
fn g0_native_fan_speed_outside_declared_domain_refuses_before_ledger_creation() {
    let fixture = fixture();
    let path = fixture.sources.join("study.fsim");
    let source = std::fs::read_to_string(&path).unwrap();
    std::fs::write(&path, source.replace(":low 0.8", ":low 0.4")).unwrap();
    let output = fs_cli::run(vec![
        "--json".into(), "study".into(),
        path.to_str().unwrap().into(),
        fixture.ledger.to_str().unwrap().into(),
    ]);
    assert_eq!(output.exit_code, fs_cli::exit::REFUSED, "{}", output.stderr);
    assert!(output.stdout.is_empty(), "admission refusal must publish no result");
    let diagnostic = JsonValue::parse(output.stderr.trim()).unwrap();
    assert_eq!(diagnostic.str_field("schema"), Some("frankensim.cli.diagnostic.v1"));
    assert_eq!(diagnostic.str_field("code"), Some("cli-uncertainty-model"));
    assert!(diagnostic.str_field("message").unwrap().contains("unchanged declared domain"));
    assert!(!fixture.ledger.exists(), "invalid fan uncertainty must fail before any solve or ledger write");
}
