//! Global effects through actual native imports, thermal solves and recovery.
use super::*;

fn source() -> String {
    STUDY.replace(":samples 4", ":samples 8")
        .replace(":method monte-carlo", ":method sobol-sensitivity")
}
fn fixture() -> Fixture {
    let f = Fixture::new();
    std::fs::write(f.sources.join("study.fsim"), source()).unwrap();
    f
}
fn report(f: &Fixture, result: &JsonValue) -> (Vec<u8>, JsonValue) {
    let receipt = result.get("receipt").unwrap();
    assert_eq!(receipt.f64_field("samples_planned"), Some(8.0));
    let bytes = artifact(&f.ledger, receipt.str_field("report_json").unwrap());
    let value = JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(value.str_field("method"), Some("sobol-sensitivity"));
    assert_eq!(value.get("statistics"), Some(&JsonValue::Null));
    for key in ["qmc", "compliance", "mean_control"] {
        assert!(value.get(key).is_none(), "pick-freeze must not publish {key}");
    }
    (bytes, value)
}
fn point(row: &JsonValue) -> Vec<f64> {
    row.get("parameters").unwrap().as_array().unwrap().iter()
        .map(|p| p.as_f64().unwrap()).collect()
}

#[test]
fn native_sensitivity_matches_hand_bound_physics_and_jansen_arithmetic() {
    let f = fixture();
    let completed = f.study(None, fs_cli::exit::SUCCESS);
    let (_, value) = report(&f, &completed);
    assert_eq!(value.str_field("termination"), Some("fixed-sensitivity-design"));
    assert_eq!(value.f64_field("evaluations_attempted"), Some(8.0));
    assert_eq!(rows(&value).len(), 8);
    let sensitivity = value.get("sensitivity").unwrap();
    assert_eq!(sensitivity.f64_field("base_rows_planned"), Some(2.0));
    assert_eq!(sensitivity.f64_field("completed_rows"), Some(2.0));
    assert_eq!(sensitivity.f64_field("evaluations_per_row"), Some(4.0));
    assert_eq!(sensitivity.f64_field("partial_row_evaluations"), Some(0.0));
    let mut outputs = Vec::new();
    for (ordinal, row) in rows(&value).iter().enumerate() {
        let (project, run, temperature) = f.independent_solve(ordinal, &point(row));
        assert_eq!(row.str_field("project_hash"), Some(project.as_str()));
        assert_eq!(row.str_field("run"), Some(run.as_str()));
        assert_eq!(row.f64_field("value_k").unwrap().to_bits(), temperature.to_bits());
        let retained = json_artifact(&f.ledger, row.str_field("qoi_receipt").unwrap());
        assert_eq!(retained.str_field("run"), Some(run.as_str()));
        outputs.push(temperature);
    }
    for block in rows(&value).chunks_exact(4) {
        let a = point(&block[0]); let b = point(&block[1]);
        for i in 0..2 {
            let mut hybrid = a.clone(); hybrid[i] = b[i];
            assert_eq!(point(&block[i + 2]), hybrid, "replace only the named physical coordinate");
        }
    }
    // An independent estimator over actual physical outputs, not a second
    // invocation of SobolExecution or statistics from the dependent hybrids.
    let bases: Vec<_> = outputs.chunks_exact(4).flat_map(|r| [r[0], r[1]]).collect();
    let mean = bases.iter().sum::<f64>() / bases.len() as f64;
    let variance = bases.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
        / (bases.len() - 1) as f64;
    assert!(variance > 0.0);
    let estimate = sensitivity.get("estimate").unwrap();
    assert!((estimate.f64_field("base_output_std_dev_k").unwrap() - variance.sqrt()).abs() < 1e-10);
    let effects = estimate.get("effects").unwrap().as_array().unwrap();
    assert_eq!(effects.len(), 2);
    for (i, (name, unit)) in [("power", "W"), ("ambient", "K")].into_iter().enumerate() {
        let total = outputs.chunks_exact(4).map(|r| (r[0] - r[i + 2]).powi(2)).sum::<f64>()
            / (4.0 * variance);
        let main = 1.0 - outputs.chunks_exact(4).map(|r| (r[1] - r[i + 2]).powi(2)).sum::<f64>()
            / (4.0 * variance);
        assert_eq!(effects[i].str_field("parameter"), Some(name));
        assert_eq!(effects[i].str_field("parameter_unit"), Some(unit));
        assert_eq!(effects[i].str_field("index_unit"), Some("1"));
        assert!((effects[i].f64_field("first_order").unwrap() - main).abs() < 1e-10);
        assert!((effects[i].f64_field("total_order").unwrap() - total).abs() < 1e-10);
    }
    let html = artifact(&f.ledger, completed.get("receipt").unwrap().str_field("report_html").unwrap());
    let html = std::str::from_utf8(&html).unwrap();
    assert!(html.contains("Main effect") && html.contains("Total effect"));
    assert!(!html.contains("empirical pass fraction:"));
}

#[test]
fn native_sensitivity_resumes_a_partial_hybrid_row_without_original_sources() {
    let full = fixture();
    let complete = full.study(None, fs_cli::exit::SUCCESS);
    let (expected_bytes, expected) = report(&full, &complete);
    let split = fixture();
    let empty = split.study(Some("0"), fs_cli::exit::BUDGET);
    let (_, empty_report) = report(&split, &empty);
    assert!(rows(&empty_report).is_empty());
    let first = split.resume(empty.str_field("run").unwrap(), Some("3"), fs_cli::exit::BUDGET);
    let (_, partial) = report(&split, &first);
    assert_eq!(rows(&partial), &rows(&expected)[..3]);
    let sensitivity = partial.get("sensitivity").unwrap();
    assert_eq!(sensitivity.f64_field("completed_rows"), Some(0.0));
    assert_eq!(sensitivity.f64_field("partial_row_evaluations"), Some(3.0));
    assert_eq!(sensitivity.get("estimate"), Some(&JsonValue::Null));
    std::fs::rename(&split.sources, split.dir.join("relocated-sensitivity-sources")).unwrap();
    let next = split.resume(first.str_field("run").unwrap(), Some("1"), fs_cli::exit::BUDGET);
    let (_, next_report) = report(&split, &next);
    assert_eq!(rows(&next_report), &rows(&expected)[..4]);
    let sensitivity = next_report.get("sensitivity").unwrap();
    assert_eq!(sensitivity.f64_field("completed_rows"), Some(1.0));
    assert_eq!(sensitivity.get("estimate"), Some(&JsonValue::Null), "one complete row is not the fixed design");
    let finished = split.resume(next.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual_bytes, _) = report(&split, &finished);
    assert_eq!(actual_bytes, expected_bytes);
    let run = finished.str_field("run").unwrap();
    let exported = command(&["--json", "report", run, split.ledger.to_str().unwrap()], fs_cli::exit::SUCCESS);
    assert_eq!(std::fs::read(exported.str_field("report_json").unwrap()).unwrap(), expected_bytes);
    assert_eq!(split.resume(run, Some("0"), fs_cli::exit::SUCCESS), finished);
}

#[test]
fn native_sensitivity_rejects_invalid_designs_before_work_and_keeps_physics_refusal_terminal() {
    let f = fixture();
    std::fs::write(f.sources.join("study.fsim"), source().replace(":samples 8", ":samples 7")).unwrap();
    f.study(None, fs_cli::exit::REFUSED);
    assert!(!f.ledger.exists(), "incomplete row budget must be refused before ledger/physics");
    std::fs::write(f.sources.join("study.fsim"), source().replace(
        ":materials (\"aa6061.fsmcdpk\")", ":materials ()")).unwrap();
    let refused = f.study(None, fs_cli::exit::REFUSED);
    let (_, value) = report(&f, &refused);
    assert!(rows(&value).is_empty());
    assert_eq!(value.f64_field("evaluations_attempted"), Some(1.0));
    assert_eq!(value.get("sensitivity").unwrap().get("estimate"), Some(&JsonValue::Null));
    assert!(value.str_field("failure").unwrap().contains("project-material-card-unknown"));
    assert_eq!(f.resume(refused.str_field("run").unwrap(), None, fs_cli::exit::REFUSED), refused);
}
