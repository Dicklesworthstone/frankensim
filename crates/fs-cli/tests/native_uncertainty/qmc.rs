use super::*;

fn source() -> String {
    STUDY.replace(":method monte-carlo",
        ":method quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
}

fn fixture() -> Fixture {
    let fixture = Fixture::new();
    std::fs::write(fixture.sources.join("study.fsim"), source()).unwrap();
    fixture
}

fn qmc(report: &JsonValue) -> &JsonValue {
    assert_eq!(report.str_field("method"), Some("quasi-monte-carlo"));
    assert_eq!(report.get("statistics"), Some(&JsonValue::Null),
        "dependent net points must not enter the iid Monte Carlo statistics");
    assert!(report.get("compliance").is_none(),
        "QMC pass proportions do not create a Bernoulli-iid confidence sequence");
    report.get("qmc").unwrap()
}

fn numbers(value: &JsonValue) -> Vec<f64> {
    value.as_array().unwrap().iter().map(|value| value.as_f64().unwrap()).collect()
}

fn between_replicate_error(means: &[f64]) -> f64 {
    let mean = means.iter().sum::<f64>() / means.len() as f64;
    (means.iter().map(|value| (value - mean).powi(2)).sum::<f64>()
        / (means.len() * (means.len() - 1)) as f64).sqrt()
}

#[test]
fn g1_native_qmc_matches_scrambled_ordinals_and_independent_physical_solves() {
    let fixture = fixture();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    let quadrature = qmc(&report);
    assert_eq!(quadrature.str_field("sampler"), Some("owen-scrambled-sobol"));
    assert_eq!(quadrature.f64_field("completed_replicates"), Some(2.0));
    assert_eq!(quadrature.f64_field("partial_replicate_samples"), Some(0.0));
    assert_eq!(quadrature.f64_field("samples_in_estimate"), Some(4.0));

    let mut plan = fs_uq::UqPlan::new("temperature-max", fs_uq::PropagationMethod::QuasiMonteCarlo, 4)
        .with_correlation(fs_uq::CorrelationModel::Independent)
        .with_compliance_threshold(297.0)
        .with_parameter(fs_uq::ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(fs_uq::ParameterUncertainty::uniform("ambient", 294.0, 300.0, "K"));
    plan.seed = 29;
    let mut sampler = fs_uq::QmcExecution::new(&plan,
        fs_uq::QmcConfig { replicates: 2, samples_per_replicate: 2 }).unwrap();
    let mut addressed = Vec::new();
    sampler.advance(4, || false, |parameters| {
        addressed.push(parameters.to_vec());
        Ok::<_, &str>(0.0)
    });
    let mut values = Vec::new();
    for (ordinal, row) in rows(&report).iter().enumerate() {
        let parameters = numbers(row.get("parameters").unwrap());
        assert_eq!(parameters.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            addressed[ordinal].iter().map(|value| value.to_bits()).collect::<Vec<_>>());
        let (project, run, value) = fixture.independent_solve(ordinal, &parameters);
        assert_eq!(row.str_field("project_hash"), Some(project.as_str()));
        assert_eq!(row.str_field("run"), Some(run.as_str()));
        assert_eq!(row.f64_field("value_k").unwrap().to_bits(), value.to_bits());
        let receipt = json_artifact(&fixture.ledger, row.str_field("qoi_receipt").unwrap());
        assert_eq!(receipt.str_field("run"), Some(run.as_str()));
        values.push(value);
    }
    assert!(values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - values.iter().copied().fold(f64::INFINITY, f64::min) > 0.1,
        "native uncertain inputs must change the solved temperature");
    let means = values.chunks_exact(2).map(|block| (block[0] + block[1]) / 2.0)
        .collect::<Vec<_>>();
    for (actual, expected) in numbers(quadrature.get("replicate_means_k").unwrap()).iter().zip(&means) {
        assert!((actual - expected).abs() < 1.0e-11);
    }
    let estimate = quadrature.get("estimate").unwrap();
    assert!((estimate.f64_field("mean_k").unwrap() - (means[0] + means[1]) / 2.0).abs() < 1.0e-11);
    assert!((estimate.f64_field("sampling_standard_error_k").unwrap()
        - between_replicate_error(&means)).abs() < 1.0e-11);
    let proportions = values.chunks_exact(2).map(|block|
        block.iter().filter(|&&value| value <= 297.0).count() as f64 / 2.0).collect::<Vec<_>>();
    let compliance = quadrature.get("compliance_estimate").unwrap();
    assert_eq!(compliance.f64_field("probability_of_compliance"),
        Some((proportions[0] + proportions[1]) / 2.0));
    assert!((compliance.f64_field("sampling_standard_error").unwrap()
        - between_replicate_error(&proportions)).abs() < 1.0e-12);
    let receipt = result.get("receipt").unwrap();
    let package = artifact(&fixture.ledger, receipt.str_field("package").unwrap());
    let text = std::str::from_utf8(&package).unwrap();
    assert!(text.contains("fixed-count-native-randomized-sobol-between-replicate-standard-error"));
    assert!(!text.contains("fixed-count-native-monte-carlo"));
    let package = fs_package::EvidencePackage::from_json(text).unwrap();
    assert!(fs_checker::check(&package).passed());
}

#[test]
fn g4_g5_native_qmc_resumes_inside_a_net_and_excludes_unfinished_points() {
    let uninterrupted = fixture();
    let complete = uninterrupted.study(None, fs_cli::exit::SUCCESS);
    let (expected_bytes, expected) = uninterrupted.report(&complete);
    let resumed = fixture();
    let zero = resumed.study(Some("0"), fs_cli::exit::BUDGET);
    let (_, empty) = resumed.report(&zero);
    assert_eq!(qmc(&empty).get("estimate"), Some(&JsonValue::Null));
    std::fs::rename(&resumed.sources, resumed.dir.join("relocated-sources")).unwrap();
    let first = resumed.resume(zero.str_field("run").unwrap(), Some("1"), fs_cli::exit::BUDGET);
    let (_, partial) = resumed.report(&first);
    assert_eq!(rows(&partial), &rows(&expected)[..1]);
    let quadrature = qmc(&partial);
    assert_eq!(quadrature.f64_field("completed_replicates"), Some(0.0));
    assert_eq!(quadrature.f64_field("partial_replicate_samples"), Some(1.0));
    assert_eq!(quadrature.f64_field("samples_in_estimate"), Some(0.0));
    assert_eq!(quadrature.get("estimate"), Some(&JsonValue::Null));
    assert_eq!(quadrature.get("compliance_estimate"), Some(&JsonValue::Null));
    let net = resumed.resume(first.str_field("run").unwrap(), Some("1"), fs_cli::exit::BUDGET);
    let (_, one_net) = resumed.report(&net);
    assert_eq!(rows(&one_net), &rows(&expected)[..2]);
    let quadrature = qmc(&one_net);
    assert_eq!(quadrature.f64_field("completed_replicates"), Some(1.0));
    assert_eq!(quadrature.f64_field("samples_in_estimate"), Some(2.0));
    assert!(quadrature.get("estimate").unwrap().f64_field("mean_k").is_some());
    assert_eq!(quadrature.get("estimate").unwrap().get("sampling_standard_error_k"),
        Some(&JsonValue::Null));
    let done = resumed.resume(net.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual_bytes, _) = resumed.report(&done);
    assert_eq!(actual_bytes, expected_bytes, "resume preserves the exact sampled native result bytes");
    for key in ["observations", "report_json", "report_html", "package", "checkpoint"] {
        assert_eq!(done.get("receipt").unwrap().str_field(key),
            complete.get("receipt").unwrap().str_field(key), "retained {key} differs");
    }
    assert_eq!(resumed.resume(done.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS), done);
}

#[test]
fn g1_native_qmc_zero_width_laws_keep_exact_solves_and_zero_replicate_variation() {
    let fixture = fixture();
    std::fs::write(fixture.sources.join("study.fsim"),
        source().replace(":low 4W :high 6W", ":low 5W :high 5W")
            .replace(":low 294K :high 300K", ":low 297K :high 297K")).unwrap();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = fixture.report(&result);
    let first = &rows(&report)[0];
    for row in rows(&report) {
        assert_eq!(numbers(row.get("parameters").unwrap()), [5.0, 297.0]);
        assert_eq!(row.str_field("run"), first.str_field("run"));
        assert_eq!(row.get("value_k"), first.get("value_k"));
    }
    let estimate = qmc(&report).get("estimate").unwrap();
    assert_eq!(estimate.get("mean_k"), first.get("value_k"));
    assert_eq!(estimate.f64_field("sampling_standard_error_k"), Some(0.0));
    assert!(report.str_field("no_claim").unwrap().contains("does not prove exactness"));
}

#[test]
fn g0_native_qmc_admission_refuses_cs_and_bad_layout_before_ledger_creation() {
    for study in [
        source().replace(":version 1", ":version 2 :compliance (bernoulli-mixture :required-probability 0.9 :alpha 0.05 :min-samples 2)"),
        source().replace(":samples-per-replicate 2", ":samples-per-replicate 3"),
        source().replace(":method quasi-monte-carlo", ":method monte-carlo"),
    ] {
        let fixture = fixture();
        std::fs::write(fixture.sources.join("study.fsim"), study).unwrap();
        let output = fs_cli::run(vec!["--json".into(), "study".into(),
            fixture.sources.join("study.fsim").to_str().unwrap().into(),
            fixture.ledger.to_str().unwrap().into()]);
        assert_eq!(output.exit_code, fs_cli::exit::REFUSED, "{}", output.stderr);
        assert!(output.stdout.is_empty());
        let diagnostic = JsonValue::parse(output.stderr.trim()).unwrap();
        assert_eq!(diagnostic.str_field("schema"), Some("frankensim.cli.diagnostic.v1"));
        assert_eq!(diagnostic.str_field("code"), Some("cli-uncertainty-model"));
        assert!(!fixture.ledger.exists(), "invalid quadrature intent must refuse before work");
    }
}
