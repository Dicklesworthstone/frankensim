//! The authored peak mode must drive the same real mechanics on run and resume.
use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-peak-stress.fsim"));
const LEGACY: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-stress.fsim"));

fn baseline(source: &str) -> Computation {
    let spec = spec::parse(source).unwrap();
    compute_observed(&spec, &CancelGate::new_clock_free(), 0, None, 0.0, |_, _| Ok(())).unwrap()
}
fn bounded(e: &StressEvaluation3) {
    assert!(e.aggregate >= e.sampled_relaxed_max);
    assert!(e.aggregate >= e.sampled_physical_max);
    assert!(e.case_relaxed_max.iter().chain(&e.case_physical_max).all(|s| *s <= e.aggregate));
    assert!(e.point_count > 0 && e.aggregate > 0.0);
}
fn calibrated_source() -> String {
    // Select a feasible, nontrivial material-removal problem from the actual
    // independently solved fixture, not a hard-coded anticipated stress value.
    let run = baseline(FIXTURE);
    assert!(run.state.audit.passed, "{:?}", run.state.error);
    let cap = 1.5 * run.state.accepted.as_ref().unwrap().aggregate;
    FIXTURE.replace(":stress-limit-pa 16.0", &format!(":stress-limit-pa {cap:.17e}"))
}
fn bytes(ledger: &Ledger, result: &Outcome, key: &str, kind: &str) -> Vec<u8> {
    let receipt = JsonValue::parse(&result.receipt).unwrap();
    linked(ledger, &receipt, key, kind).unwrap()
}
fn document(bytes: &[u8]) -> JsonValue {
    JsonValue::parse(std::str::from_utf8(bytes).unwrap()).unwrap()
}

#[test]
fn peak_type_is_explicit_identity_bearing_and_zero_weights_do_not_disable_cases() {
    let spec = spec::parse(FIXTURE).unwrap();
    assert_eq!(spec.stress_measure, StressMeasure3::SampledPeakBound);
    let again = spec::parse(&spec.canonical).unwrap();
    assert_eq!(again.id, spec.id);
    assert_eq!(again.stress_measure, spec.stress_measure);
    let average = spec::parse(&FIXTURE.replace("sampled-peak-limited-simp", "stress-limited-simp")).unwrap();
    assert_eq!(average.stress_measure, StressMeasure3::NormalizedAverage);
    assert_ne!(average.id, spec.id);
    assert_eq!(spec::parse(LEGACY).unwrap().stress_measure, StressMeasure3::NormalizedAverage);
    assert!(spec::parse(&FIXTURE.replace(":weight 0.3", ":weight 0.0")).is_ok());
    assert!(spec::parse(&LEGACY.replace(":weight 0.3", ":weight 0.0")).is_err());
    for source in [
        FIXTURE.replace(":weight 0.3", ":weight 0.0").replace(":weight 0.7", ":weight 0.0"),
        FIXTURE.replace(":weight 0.3", ":weight -1.0"),
        FIXTURE.replace(":weight 0.3", ":weight 1.1"),
        FIXTURE.replace("sampled-peak-limited-simp", "continuum-peak-limited-simp"),
        FIXTURE.replace(":aggregation-power 16.0", ":aggregation-power 65.0"),
        FIXTURE.replace(":density-floor 0.05", ":density-floor 0.0001"),
        FIXTURE.replace(":beta 0.0", ":beta 0.0 :ignore-worst-case true"),
    ] { assert!(spec::parse(&source).is_err(), "must not ignore an unsupported declaration"); }
}

#[test]
fn hot_low_weight_native_case_cannot_be_hidden_by_the_normalized_average() {
    let hot = FIXTURE.replace("(0.0 0.0 -1.0)", "(0.0 0.0 -20.0)")
        .replace(":weight 0.7", ":weight 1e-30");
    let average = baseline(&hot.replace("sampled-peak-limited-simp", "stress-limited-simp"));
    let peak = baseline(&hot);
    let zero = baseline(&hot.replace(":weight 1e-30", ":weight 0.0"));
    let a = average.state.accepted.as_ref().unwrap();
    let b = peak.state.accepted.as_ref().unwrap();
    let c = zero.state.accepted.as_ref().unwrap();
    assert!(peak.state.audit.passed && zero.state.audit.passed);
    bounded(b); bounded(c);
    assert!(a.case_physical_max[1] > a.aggregate);
    assert!(b.aggregate > 5.0 * a.aggregate);
    assert_eq!(b.aggregate.to_bits(), c.aggregate.to_bits());
    assert_eq!(b.gradient, c.gradient);
    assert_eq!(b.displacements, c.displacements);
    assert_eq!(b.displacements, a.displacements, "the constraint may not change the applied loads");
    assert!(c.adjoints[1].iter().any(|z| z.abs() > 0.0));
}

#[test]
fn peak_native_updates_pass_gradients_and_select_a_sample_feasible_reduced_design() {
    let source = calibrated_source();
    let spec = spec::parse(&source).unwrap();
    let mut published = Vec::new();
    let run = compute_observed(&spec, &CancelGate::new_clock_free(), 4, None, 0.0, |_, state| {
        assert!(state.audit.passed);
        bounded(state.accepted.as_ref().unwrap());
        published.push(state.iterations());
        Ok(())
    }).unwrap();
    assert!(published.len() > 1, "must publish real accepted updates");
    assert_eq!(run.state.iterations(), 4);
    assert!(run.state.audit.passed, "{:?}", run.state.error);
    let best = run.state.best.as_ref().unwrap();
    bounded(best);
    assert!(best.volume_fraction < run.state.history[0].volume_fraction);
    let options = spec.stress.unwrap();
    let accepted_cap = options.stress_limit * (1.0 + options.optimizer.tolerance);
    assert!(best.sampled_physical_max <= accepted_cap);
    assert!(best.sampled_relaxed_max <= accepted_cap);
    assert_ne!(best.displacements[0], best.displacements[1]);
}

#[test]
fn peak_ledger_resume_preserves_exact_fields_policy_and_solver_free_exports() {
    let source = calibrated_source();
    let spec = spec::parse(&source).unwrap();
    let gate = CancelGate::new_clock_free();
    let full_db = Ledger::open(":memory:").unwrap();
    let full = drive(&spec, &full_db, Some(4), &gate).unwrap();
    let split_db = Ledger::open(":memory:").unwrap();
    let first = drive(&spec, &split_db, Some(2), &gate).unwrap();
    let old = load(&split_db, &first.pointer).unwrap();
    assert_eq!(integer(&old.value, "iterations_completed").unwrap(), 2);
    let resumed = super::resume(&split_db, &old, Some(2), &gate).unwrap();
    let receipt = JsonValue::parse(&resumed.receipt).unwrap();
    assert_eq!(integer(&receipt, "iterations_completed").unwrap(), 4);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(bytes(&full_db, &full, key, kind), bytes(&split_db, &resumed, key, kind));
    }
    let summary = document(&bytes(&split_db, &resumed, "report_json", "study-report-json"));
    let measure = summary.get("stress_measure").unwrap();
    assert_eq!(measure.str_field("kind"), Some("unweighted-sampled-peak-bound"));
    assert_eq!(measure.get("sampled_maximum_stress_is_constrained"), Some(&JsonValue::Bool(true)));
    assert_eq!(measure.get("continuum_maximum_stress_is_constrained"), Some(&JsonValue::Bool(false)));
    assert_eq!(measure.get("load_weights_affect_constraint"), Some(&JsonValue::Bool(false)));
    let html = bytes(&split_db, &resumed, "report_html", "study-report-html");
    assert!(std::str::from_utf8(&html).unwrap().contains("unweighted sampled-peak bound"));
    let package = bytes(&split_db, &resumed, "package", "study-package");
    assert!(std::str::from_utf8(&package).unwrap().contains("unweighted-sampled-peak-bound"));
    let full_work = Spent::read(&JsonValue::parse(&full.receipt).unwrap()).unwrap();
    let resumed_work = Spent::read(&receipt).unwrap();
    assert!(resumed_work.linear.linear_solves > full_work.linear.linear_solves,
        "rebuilding physical endpoints is charged, not refunded");
    let changed = spec::parse(&source.replace("sampled-peak-limited-simp", "stress-limited-simp")).unwrap();
    let error = match resume::drive(&changed, &split_db, Some(1), &gate, Some(&old)) {
        Err(error) => error,
        Ok(_) => panic!("cannot authorize a different stress functional using a retained peak checkpoint"),
    };
    assert_eq!(error.code, "cli-study-sdf3-stress-checkpoint");
}

#[test]
fn native_peak_prescribed_void_retains_the_modeled_physical_stress() {
    let source = format!("{}\n(design-regions :solid () :void (((0.0 0.0 0.0) (1.0 1.0 1.0)))))",
        FIXTURE.trim_end().strip_suffix(')').unwrap());
    let run = baseline(&source);
    assert!(run.state.audit.passed, "{:?}", run.state.error);
    let e = run.state.accepted.as_ref().unwrap();
    assert!(e.projected_rho.iter().all(|r| *r == 0.0));
    assert_eq!(e.sampled_relaxed_max, 0.0);
    assert!(e.sampled_physical_max > 0.0, "ersatz material still carries the prescribed loads");
    bounded(e);
}

#[test]
fn native_peak_work_exhaustion_and_precancellation_publish_no_unchecked_design() {
    for source in [
        FIXTURE.replace(":max-evaluations 2000", ":max-evaluations 1"),
        FIXTURE.replace(":linear-iterations 250000", ":linear-iterations 1"),
        FIXTURE.replace(":max-stress-points 50000", ":max-stress-points 1"),
    ] {
        let spec = spec::parse(&source).unwrap();
        let run = compute_observed(&spec, &CancelGate::new_clock_free(), 1, None, 0.0, |_, _| {
            panic!("no publication before a complete gradient-checked baseline")
        }).unwrap();
        assert_eq!(run.state.status, "budget-exhausted");
        assert!(run.state.selected().is_none());
        assert!(run.state.spent.linear.linear_iterations <= spec.linear);
    }
    let spec = spec::parse(FIXTURE).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    let error = match compute_observed(&spec, &gate, 1, None, 0.0, |_, _| panic!("cancelled")) {
        Err(error) => error,
        Ok(_) => panic!("precancelled study cannot start geometry"),
    };
    assert_eq!(error.exit, exit::CANCELLED);
}
