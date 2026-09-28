use super::*;

fn source(qmc: bool, rho: f64) -> String {
    let declared = STUDY.replace(":correlation independent", &format!(
        ":correlation (gaussian-copula :latent-correlation ((1 {rho}) ({rho} 1)))"));
    if qmc {
        declared.replace(":method monte-carlo",
            ":method quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
    } else { declared }
}

fn fixture(source: &str) -> Fixture {
    let fixture = Fixture::new();
    std::fs::write(fixture.sources.join("study.fsim"), source).unwrap();
    fixture
}

fn verify_law(fixture: &Fixture, result: &JsonValue, report: &JsonValue, rho: f64) {
    assert_eq!(report.str_field("correlation"), Some("gaussian-copula"));
    let dependence = report.get("dependence").unwrap();
    assert_eq!(dependence.str_field("matrix_coordinates"), Some("latent-standard-normal"));
    let order = dependence.get("parameter_order").unwrap().as_array().unwrap();
    assert_eq!(order.iter().map(|value| value.as_str().unwrap()).collect::<Vec<_>>(), ["power", "ambient"]);
    let matrix = dependence.get("latent_correlation").unwrap().as_array().unwrap();
    assert_eq!(matrix[0].as_array().unwrap()[1].as_f64(), Some(rho));
    assert_eq!(matrix[1].as_array().unwrap()[0].as_f64(), Some(rho));
    assert!(report.str_field("no_claim").unwrap().contains("not physical Pearson correlation"));
    let receipt = result.get("receipt").unwrap();
    for key in ["report_html", "package"] {
        let bytes = artifact(&fixture.ledger, receipt.str_field(key).unwrap());
        let text = std::str::from_utf8(&bytes).unwrap();
        assert!(text.contains("Gaussian copula"), "{key} omitted the declared law");
        assert!(!text.contains("independent uniform"), "{key} mislabeled dependent inputs");
    }
}

#[test]
fn g1_native_copula_mc_and_qmc_map_dependent_inputs_into_real_solves() {
    for (qmc, rho) in [(false, 1.0), (true, -1.0)] {
        let fixture = fixture(&source(qmc, rho));
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = fixture.report(&result);
        verify_law(&fixture, &result, &report, rho);
        assert_eq!(report.str_field("method"), Some(if qmc { "quasi-monte-carlo" } else { "monte-carlo" }));
        assert!(report.get("compliance").is_none());
        for (ordinal, row) in rows(&report).iter().enumerate() {
            let parameters = row.get("parameters").unwrap().as_array().unwrap().iter()
                .map(|value| value.as_f64().unwrap()).collect::<Vec<_>>();
            let power_quantile = (parameters[0] - 4.0) / 2.0;
            let ambient_quantile = (parameters[1] - 294.0) / 6.0;
            assert!((0.0..=1.0).contains(&power_quantile));
            assert!((0.0..=1.0).contains(&ambient_quantile));
            let expected = if rho == 1.0 { power_quantile } else { 1.0 - power_quantile };
            assert!((ambient_quantile - expected).abs() < 1.0e-12,
                "singular copula dependence must survive the physical-unit transform");
            let (project, run, value) = fixture.independent_solve(ordinal, &parameters);
            assert_eq!(row.str_field("project_hash"), Some(project.as_str()));
            assert_eq!(row.str_field("run"), Some(run.as_str()));
            assert_eq!(row.f64_field("value_k").unwrap().to_bits(), value.to_bits());
        }
        let values = rows(&report).iter().map(|row| row.f64_field("value_k").unwrap()).collect::<Vec<_>>();
        assert!(values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - values.iter().copied().fold(f64::INFINITY, f64::min) > 0.1);
        if qmc {
            assert_eq!(report.get("statistics"), Some(&JsonValue::Null));
            assert_eq!(report.get("qmc").unwrap().f64_field("completed_replicates"), Some(2.0));
        } else {
            assert!(report.get("statistics").unwrap().f64_field("mean_k").is_some());
        }
    }
}

#[test]
fn g4_g5_native_copula_replay_retains_the_joint_law_and_physical_sample_prefix() {
    for qmc in [false, true] {
        let declared = source(qmc, 0.65);
        let full = fixture(&declared);
        let completed = full.study(None, fs_cli::exit::SUCCESS);
        let (expected_bytes, expected) = full.report(&completed);
        let resumed = fixture(&declared);
        let first = resumed.study(Some("1"), fs_cli::exit::BUDGET);
        let (_, prefix) = resumed.report(&first);
        assert_eq!(rows(&prefix), &rows(&expected)[..1]);
        if qmc {
            assert_eq!(prefix.get("qmc").unwrap().get("estimate"), Some(&JsonValue::Null));
        }
        // An edited source must not silently change a retained run's joint law.
        std::fs::write(resumed.sources.join("study.fsim"), STUDY).unwrap();
        std::fs::rename(&resumed.sources, resumed.dir.join("relocated-sources")).unwrap();
        let finished = resumed.resume(first.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
        let (actual_bytes, report) = resumed.report(&finished);
        assert_eq!(actual_bytes, expected_bytes);
        verify_law(&resumed, &finished, &report, 0.65);
        for key in ["observations", "report_json", "report_html", "package", "checkpoint"] {
            assert_eq!(finished.get("receipt").unwrap().str_field(key),
                completed.get("receipt").unwrap().str_field(key), "{key} differs after resume");
        }
    }
}

#[test]
fn g1_native_copula_bernoulli_policy_uses_raw_iid_vector_outcomes() {
    let declared = source(false, 0.8).replace(":version 1", ":version 2")
        .replace(":low 4W :high 6W", ":low 5W :high 5W")
        .replace(":parameters (", ":compliance (bernoulli-mixture :required-probability 0.5 :alpha 0.05 :min-samples 4) :parameters (");
    let fixture = fixture(&declared);
    let result = fixture.study(None, fs_cli::exit::BUDGET);
    let (_, report) = fixture.report(&result);
    assert_eq!(report.get("statistics"), Some(&JsonValue::Null));
    assert_eq!(super::compliance::verify_confidence(&report, 4, 0.5), None);
    assert!(rows(&report).iter().all(|row|
        row.get("parameters").unwrap().as_array().unwrap()[0].as_f64() == Some(5.0)),
        "a zero-width marginal remains constant without removing its latent coordinate");
    verify_law(&fixture, &result, &report, 0.8);
    assert_eq!(report.str_field("termination"), Some("lifetime-sample-budget"));
}

#[test]
fn g0_native_copula_numerical_admission_precedes_ledger_creation() {
    let indefinite = source(false, 0.8)
        .replace("((1 0.8) (0.8 1))", "((1 0.9 0.9) (0.9 1 -0.9) (0.9 -0.9 1))")
        .replace(":parameters (", ":parameters ((uniform :name \"fan\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.6 :high 1.4)");
    for declared in [
        source(false, 0.8).replace("((1 0.8) (0.8 1))", "((1 0.2) (0.3 1))"),
        source(false, 2.0),
        indefinite.clone(),
        indefinite.replace(":method monte-carlo",
            ":method quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)"),
    ] {
        let fixture = fixture(&declared);
        let output = fs_cli::run(vec!["--json".into(), "study".into(),
            fixture.sources.join("study.fsim").to_str().unwrap().into(),
            fixture.ledger.to_str().unwrap().into()]);
        assert_eq!(output.exit_code, fs_cli::exit::REFUSED, "{}", output.stderr);
        assert!(output.stdout.is_empty());
        let diagnostic = JsonValue::parse(output.stderr.trim()).unwrap();
        assert_eq!(diagnostic.str_field("code"), Some("cli-uncertainty-plan"));
        assert!(output.stderr.contains("correlation"));
        assert!(!fixture.ledger.exists(), "invalid joint law must refuse before physical work");
    }
}
