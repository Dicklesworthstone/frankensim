use super::*;

const FIXTURE: &str = include_str!("../../../../../../examples/marquee/bracket-3d-load-space.fsim");

fn compute(spec: &Spec, updates: usize) -> Computation {
    compute_observed(spec, &CancelGate::new(), updates, None, 0.0, |_, _| Ok(())).unwrap()
}

#[test]
fn g0_load_space_schema_binds_the_force_set_and_physical_cap() {
    let spec = spec::parse(FIXTURE).unwrap();
    assert!(spec.stress.is_none());
    assert!(spec.load_space.is_some());
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    assert!(ir_for(LOAD_SPACE3_DRIVER, spec.id, 0).contains("\"objective\":\"1\""));
    for source in [
        FIXTURE.replace(":weight 1.0", ":weight 0.5"),
        FIXTURE.replace(":compliance-limit-j 10.0", ":compliance-limit-j 0.0"),
        FIXTURE.replace(":spectral-temperature 0.05", ":spectral-temperature 0.0"),
        FIXTURE.replace(":threshold-offset 0.1", ":threshold-offset 0.5"),
        FIXTURE.replace(":maximum-level 1", ":maximum-level 2"),
        FIXTURE.replace(":initial-density 0.75", ":initial-density 0.01"),
        FIXTURE.replace(
            ":max-optimizer-evaluations 2000",
            ":max-optimizer-evaluations 0",
        ),
        FIXTURE.replace(":unit \"1\"", ":unit \"J\""),
        FIXTURE.replace(":type dilated-volume-fraction", ":type volume-fraction"),
    ] {
        assert!(
            spec::parse(&source).is_err(),
            "must reject an ambiguous or invalid load-space declaration"
        );
    }
    for source in [
        FIXTURE.replace(":compliance-limit-j 10.0", ":compliance-limit-j 11.0"),
        FIXTURE.replace(":threshold-offset 0.1", ":threshold-offset 0.2"),
        FIXTURE.replace("(0.0 -1.0 -0.5)", "(0.0 1.0 -0.5)"),
    ] {
        assert_ne!(spec::parse(&source).unwrap().id, spec.id);
    }
}

#[test]
fn g1_g3_simultaneous_force_combinations_can_fail_when_each_basis_load_passes() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    let initial = compute(&spec, 0);
    let e = &initial.state.accepted.evaluation;
    let basis_max = e
        .physical
        .eroded
        .objective
        .case_compliances
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    assert!(
        e.worst_compliance > 1.1 * basis_max,
        "correlated basis columns must expose cross-work"
    );
    let cap = 0.5 * (basis_max + e.worst_compliance);
    spec.load_space.as_mut().unwrap().cap_j = cap;
    let hostile = compute(&spec, 0);
    let evaluated = &hostile.state.accepted.evaluation;
    assert!(
        evaluated
            .physical
            .eroded
            .objective
            .case_compliances
            .iter()
            .all(|&c| c < cap)
    );
    assert!(evaluated.worst_compliance > cap);
    assert!(hostile.state.best.is_none());
    assert_ne!(hostile.state.status, "completed");
    let report = output::measurement(&hostile.state.accepted, spec.load_space.unwrap()).unwrap();
    assert_eq!(
        JsonValue::parse(&report).unwrap().get("feasible"),
        Some(&JsonValue::Bool(false))
    );
}

#[test]
fn g1_g5_native_load_space_design_reduces_material_with_real_adjoint_gates() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    let initial = compute(&spec, 0);
    let baseline = initial
        .state
        .accepted
        .evaluation
        .physical
        .dilated_volume_fraction;
    spec.load_space.as_mut().unwrap().cap_j =
        1.5 * initial.state.accepted.evaluation.worst_compliance;
    let run = compute(&spec, 8);
    assert!(run.state.gradient.numerical.passed());
    assert_eq!(run.state.gradient.numerical.probes.len(), 2);
    let best = run
        .state
        .best
        .as_ref()
        .expect("retain the least-material feasible accepted point");
    let config = spec.load_space.unwrap();
    assert!(best.evaluation.physical.dilated_volume_fraction < baseline);
    assert!(
        config
            .cap(spec.loads.len())
            .unwrap()
            .admits(&best.evaluation, config.options.optimizer.tolerance)
    );
    let e = &best.evaluation.physical;
    assert!(
        e.eroded
            .projected_rho
            .iter()
            .zip(&e.nominal_rho)
            .zip(&e.dilated_rho)
            .all(|((&a, &b), &c)| a <= b + 1e-14 && b <= c + 1e-14)
    );
    assert!(run.state.spent.linear.linear_iterations > 0);
    assert!(run.state.optimizer.evaluations <= config.options.optimizer.max_evaluations);
    let replay = compute(&spec, 8);
    assert_eq!(replay.state.history, run.state.history);
    assert_eq!(replay.state.selected().rho, run.state.selected().rho);
    assert_eq!(
        replay.state.selected().evaluation.worst_displacement,
        best.evaluation.worst_displacement
    );
    assert_eq!(replay.state.spent.linear, run.state.spent.linear);
}

#[test]
fn g4_native_load_space_stops_respect_cancellation_and_original_work_caps() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    let gate = CancelGate::new();
    gate.request();
    let error = compute_observed(&spec, &gate, 1, None, 0.0, |_, _| {
        panic!("no accepted output before geometry")
    })
    .err()
    .expect("precancellation must refuse");
    assert_eq!(error.exit, exit::CANCELLED);
    spec.linear = 1;
    let error = compute_observed(&spec, &CancelGate::new(), 1, None, 0.0, |_, _| {
        panic!("no admitted baseline")
    })
    .err()
    .expect("linear allowance cannot admit the gradient gate");
    assert_eq!(error.exit, exit::BUDGET);
    let mut spec = spec::parse(FIXTURE).unwrap();
    spec.load_space
        .as_mut()
        .unwrap()
        .options
        .optimizer
        .max_evaluations = 1;
    let run = compute(&spec, 8);
    assert_eq!(run.state.optimizer.evaluations, 1);
    assert_eq!(run.state.optimizer.iterations, 0);
    assert_eq!(run.state.status, "budget-exhausted");
    assert_eq!(run.state.history.len(), 1);

    let spec = spec::parse(FIXTURE).unwrap();
    let gate = CancelGate::new();
    let mut retained = 0;
    let stopped = compute_observed(&spec, &gate, 4, None, 0.0, |study, state| {
        assert_eq!(study.cells(), state.accepted.rho.len());
        assert!(state.gradient.numerical.passed());
        retained += 1;
        gate.request();
        Ok(())
    })
    .unwrap();
    assert_eq!(retained, 1);
    assert_eq!(stopped.state.status, "cancelled");
    assert_eq!(stopped.state.optimizer.iterations, 0);
    assert_eq!(stopped.state.history.len(), 1);
    assert!(stopped.state.spent.linear.linear_iterations > 0);
}

#[test]
fn g3_native_load_space_retains_real_designs_for_standard_report_and_package_exports() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, Some(1), &CancelGate::new()).unwrap();
    let retained = load(&ledger, &out.pointer).unwrap();
    assert_eq!(retained.value.str_field("driver"), Some(LOAD_SPACE3_DRIVER));
    assert_eq!(
        retained.value.get("resume_supported"),
        Some(&JsonValue::Bool(true))
    );
    let summary = linked(&ledger, &retained.value, "report_json", "study-report-json").unwrap();
    let report = JsonValue::parse(std::str::from_utf8(&summary).unwrap()).unwrap();
    assert_eq!(
        report.str_field("objective"),
        Some("dilated-volume-fraction")
    );
    assert!(
        report
            .get("selected")
            .unwrap()
            .get("worst_coefficients")
            .unwrap()
            .as_array()
            .is_some()
    );
    let design = linked(&ledger, &retained.value, "design", "study-design").unwrap();
    let design = JsonValue::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    assert!(
        design
            .get("selected")
            .unwrap()
            .get("basis_displacements_m")
            .unwrap()
            .as_array()
            .is_some()
    );
    let package = linked(&ledger, &retained.value, "package", "study-package").unwrap();
    let package = EvidencePackage::from_json(std::str::from_utf8(&package).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
    let continued = super::super::resume(&ledger, &retained, Some("1"), &CancelGate::new()).unwrap();
    let continued = load(&ledger, &continued.pointer).unwrap();
    assert!(integer(&continued.value, "iterations_completed").unwrap()
        >= integer(&retained.value, "iterations_completed").unwrap());
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, retained.bytes);
}
