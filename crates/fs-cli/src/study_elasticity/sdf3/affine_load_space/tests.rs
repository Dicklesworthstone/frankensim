use super::*;

const FIXTURE: &str =
    include_str!("../../../../../../examples/marquee/bracket-3d-affine-load-space.fsim");

fn fixture() -> Spec {
    spec::parse(FIXTURE).unwrap()
}

fn compute(spec: &Spec, updates: usize) -> Computation {
    compute_observed(spec, &CancelGate::new(), updates, None, 0.0, |_, _| Ok(())).unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1e-20),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_affine_witness_keeps_self_weight_out_of_the_uncertain_columns() {
    let spec = fixture();
    let run = compute(&spec, 0);
    let e = &run.state.accepted.evaluation;
    assert!(run.state.gradient.passed());
    assert!(e.oracle.gradient_admissible());
    assert!(e.offset_force.iter().any(|f| *f != 0.0));
    assert!(e.linear.iter().any(|f| *f != 0.0));
    let law = spec.material_acceleration.unwrap();
    let op = run.study.operator().elasticity();
    let body = op
        .prepare_material_body_load(&|_| law.acceleration_m_s2, law.max_terms, || {
            ControlFlow::Continue(())
        })
        .unwrap();
    let offset = body
        .load(&e.physical.eroded_rho, || ControlFlow::Continue(()))
        .unwrap();
    for (actual, expected) in e.offset_force.iter().zip(&offset) {
        near(*actual, *expected);
    }
    loading::with_laws(&spec, |laws| {
        let basis: Vec<_> = laws
            .iter()
            .map(|law| {
                op.reference_load(law.load, || ControlFlow::Continue(()))
                    .unwrap()
            })
            .collect();
        for d in 0..e.worst_force.len() {
            let f = offset[d]
                + basis
                    .iter()
                    .zip(&e.worst_coefficients)
                    .map(|(f, q)| f[d] * q)
                    .sum::<f64>();
            let u = e.offset_displacement[d]
                + e.basis_displacements
                    .iter()
                    .zip(&e.worst_coefficients)
                    .map(|(u, q)| u[d] * q)
                    .sum::<f64>();
            near(e.worst_force[d], f);
            near(e.worst_displacement[d], u);
        }
    });
    near(
        e.worst_compliance,
        e.worst_force
            .iter()
            .zip(&e.worst_displacement)
            .map(|(f, u)| f * u)
            .sum(),
    );
    assert!(e.witness_relative_residual < 1e-7);
    assert!(e.oracle.relative_duality_gap < 1e-8);
}

#[test]
fn g1_g5_affine_design_reduces_material_under_the_actual_total_load_cap() {
    let mut spec = fixture();
    let initial = compute(&spec, 0);
    let baseline = initial
        .state
        .accepted
        .evaluation
        .physical
        .dilated_volume_fraction;
    spec.affine_load_space.as_mut().unwrap().cap_j =
        1.5 * initial.state.accepted.evaluation.worst_compliance;
    let run = compute(&spec, 8);
    let best = run
        .state
        .best
        .as_ref()
        .expect("initial feasible point remains available");
    assert!(best.evaluation.physical.dilated_volume_fraction < baseline);
    assert!(spec.affine_load_space.unwrap().feasible(best));
    assert!(
        run.state.optimizer.evaluations
            <= spec
                .affine_load_space
                .unwrap()
                .options
                .optimizer
                .max_evaluations
    );
    assert!(run.state.spent.linear.linear_iterations > 0);
    assert!(run.state.spent.linear.linear_iterations <= spec.linear);
    let replay = compute(&spec, 8);
    assert_eq!(run.state.history, replay.state.history);
    assert_eq!(run.state.selected().rho, replay.state.selected().rho);
    assert_eq!(
        run.state.selected().evaluation.worst_displacement,
        replay.state.selected().evaluation.worst_displacement
    );
    let mut infeasible = spec.clone();
    infeasible.affine_load_space.as_mut().unwrap().cap_j =
        0.5 * initial.state.accepted.evaluation.worst_compliance;
    let refused_cap = compute(&infeasible, 0);
    assert!(refused_cap.state.best.is_none());
    assert!(
        !infeasible
            .affine_load_space
            .unwrap()
            .feasible(&refused_cap.state.accepted)
    );
    assert_ne!(refused_cap.state.status, "completed");
}

#[test]
fn g4_affine_design_stops_without_losing_accepted_fields_or_renewing_work() {
    let mut spec = fixture();
    let gate = CancelGate::new();
    gate.request();
    let stopped = compute_observed(&spec, &gate, 1, None, 0.0, |_, _| {
        panic!("pre-cancelled work cannot publish a baseline")
    });
    assert_eq!(stopped.err().unwrap().exit, exit::CANCELLED);
    spec.linear = 1;
    let stopped = compute_observed(&spec, &CancelGate::new(), 1, None, 0.0, |_, _| {
        panic!("one Krylov iteration cannot pass the gradient gate")
    });
    assert_eq!(stopped.err().unwrap().exit, exit::BUDGET);
    let mut spec = fixture();
    spec.affine_load_space
        .as_mut()
        .unwrap()
        .options
        .optimizer
        .max_evaluations = 1;
    let run = compute(&spec, 8);
    assert_eq!(run.state.optimizer.evaluations, 1);
    assert_eq!(run.state.optimizer.iterations, 0);
    assert_eq!(run.state.status, "budget-exhausted");
    let spec = fixture();
    let gate = CancelGate::new();
    let mut retained = None;
    let run = compute_observed(&spec, &gate, 4, None, 0.0, |_, state| {
        retained = Some(state.accepted.rho.clone());
        gate.request();
        Ok(())
    })
    .unwrap();
    assert_eq!(run.state.status, "cancelled");
    assert_eq!(Some(run.state.accepted.rho), retained);
    assert_eq!(run.state.optimizer.iterations, 0);
    assert!(run.state.spent.linear.linear_iterations > 0);
}

#[test]
fn g3_affine_artifacts_keep_the_physical_offset_and_optimizer_checkpoint() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, Some(0), &CancelGate::new()).unwrap();
    let retained = load(&ledger, &out.pointer).unwrap();
    assert_eq!(
        retained.value.str_field("driver"),
        Some(AFFINE_LOAD_SPACE3_DRIVER)
    );
    assert_eq!(
        retained.value.get("resume_supported"),
        Some(&JsonValue::Bool(true))
    );
    let report = linked(&ledger, &retained.value, "report_json", "study-report-json").unwrap();
    let report = JsonValue::parse(std::str::from_utf8(&report).unwrap()).unwrap();
    assert!(
        report
            .get("selected")
            .unwrap()
            .get("normalized_offset")
            .unwrap()
            .as_f64()
            .unwrap()
            > 0.0
    );
    let design = linked(&ledger, &retained.value, "design", "study-design").unwrap();
    let fields = JsonValue::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    assert!(
        fields
            .get("selected")
            .unwrap()
            .get("offset_displacement_m")
            .unwrap()
            .as_array()
            .is_some()
    );
    let package = linked(&ledger, &retained.value, "package", "study-package").unwrap();
    let package = EvidencePackage::from_json(std::str::from_utf8(&package).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
    assert_eq!(retained.value.str_field("resume_mode"), Some(resume::MODE));
    let checkpoint = linked(&ledger, &retained.value, "checkpoint", resume::KIND).unwrap();
    assert!(!checkpoint.is_empty());
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, retained.bytes);
}
