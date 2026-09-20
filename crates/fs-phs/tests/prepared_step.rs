use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage, step};

fn oscillator(damping: f64) -> PortHamiltonian {
    PortHamiltonian::new(2, 1, vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping], vec![0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![9.0, 0.0, 0.0, 1.0], 2).unwrap())).unwrap()
}

#[derive(Debug)]
struct CoupledQuartic;
impl Storage for CoupledQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let stretch = x[0] * x[0] + 0.3 * x[2] * x[2];
        0.5 * (7.0 * x[0] * x[0] + x[1] * x[1] + 11.0 * x[2] * x[2] + x[3] * x[3])
            + 20.0 * stretch * stretch
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let stretch = x[0] * x[0] + 0.3 * x[2] * x[2];
        out[0] = 7.0 * x[0] + 80.0 * stretch * x[0];
        out[1] = x[1];
        out[2] = 11.0 * x[2] + 24.0 * stretch * x[2];
        out[3] = x[3];
    }
}
fn nonlinear() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0; j[4] = -1.0; j[11] = 1.0; j[14] = -1.0;
    PortHamiltonian::new(4, 1, j, vec![0.0; 16], vec![0.0, 1.0, 0.0, 0.4], Box::new(CoupledQuartic)).unwrap()
}

#[test]
fn prepared_step_matches_analytic_implicit_midpoint() {
    let sys = oscillator(0.0);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut x, mut y) = ([0.0; 2], [0.0]);
    let dt = 0.01;
    let record = work.step_into(&sys, &[0.3, 0.7], &[0.0], dt, &mut x, &mut y).unwrap();
    let denominator = 1.0 + 0.25 * dt * dt * 9.0;
    let q = ((1.0 - 0.25 * dt * dt * 9.0) * 0.3 + dt * 0.7) / denominator;
    let p = (-dt * 9.0 * 0.3 + (1.0 - 0.25 * dt * dt * 9.0) * 0.7) / denominator;
    assert!((x[0] - q).abs() < 1e-11);
    assert!((x[1] - p).abs() < 1e-11);
    assert!((y[0] - 0.5 * (0.7 + p)).abs() < 1e-10);
    assert!(record.balance_residual().abs() < 1e-11);
}

#[test]
fn prepared_step_retains_nonlinear_mode_exchange_and_energy() {
    let sys = nonlinear();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut actual = [0.25, 0.4, -0.13, 0.1];
    let mut reference = actual.to_vec();
    let initial_energy = sys.hamiltonian(&actual);
    let mut next = [0.0; 4];
    let mut y = [0.0];
    for _ in 0..400 {
        let expected = step(&sys, &reference, &[0.0], 0.002).unwrap();
        let record = work.step_into(&sys, &actual, &[0.0], 0.002, &mut next, &mut y).unwrap();
        for (a, b) in next.iter().zip(&expected.x) { assert!((a - b).abs() < 2e-8); }
        assert!((y[0] - expected.y[0]).abs() < 2e-8);
        assert!(record.balance_residual().abs() < 1e-9);
        assert!((sys.hamiltonian(&next) - initial_energy).abs() < 1e-8);
        actual = next;
        reference = expected.x;
    }
}

#[test]
fn prepared_step_accounts_for_external_work_and_dissipation() {
    let sys = oscillator(0.7);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut state = [0.2, -0.4];
    let mut next = [0.0; 2]; let mut y = [0.0];
    for i in 0..100 {
        let force = if i < 30 { 0.8 } else { 0.0 };
        let record = work.step_into(&sys, &state, &[force], 0.003, &mut next, &mut y).unwrap();
        assert!((record.supplied - 0.003 * force * y[0]).abs() < 1e-14);
        assert!(record.dissipated >= 0.0);
        assert!(record.supply_defect() <= 1e-10);
        assert!(record.balance_residual().abs() < 1e-10);
        state = next;
    }
}

#[test]
fn prepared_step_failures_preserve_both_outputs_and_do_not_poison_retry() {
    let sys = oscillator(0.0);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([123.0, 456.0], [789.0]);
    for (state, input, dt) in [([f64::NAN, 0.0], [0.0], 0.01), ([0.0, 0.0], [f64::INFINITY], 0.01), ([0.0, 0.0], [0.0], -0.01)] {
        assert!(work.step_into(&sys, &state, &input, dt, &mut next, &mut y).is_err());
        assert_eq!(next, [123.0, 456.0]); assert_eq!(y, [789.0]);
    }
    assert!(work.step_into(&sys, &[0.0], &[0.0], 0.01, &mut next, &mut y).is_err());
    work.set_iteration_limit(0).unwrap();
    assert!(work.step_into(&sys, &[0.3, 0.7], &[0.0], 0.01, &mut next, &mut y).is_err());
    assert_eq!(next, [123.0, 456.0]); assert_eq!(y, [789.0]);
    work.set_iteration_limit(50).unwrap();
    work.step_into(&sys, &[0.3, 0.7], &[0.0], 0.01, &mut next, &mut y).unwrap();
    let expected = step(&sys, &[0.3, 0.7], &[0.0], 0.01).unwrap();
    for (a, b) in next.iter().zip(expected.x) { assert!((a - b).abs() < 1e-10); }
}

#[test]
fn zero_duration_is_a_noop_with_the_correct_port_output() {
    let sys = oscillator(0.5);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([0.0; 2], [0.0]);
    let record = work.step_into(&sys, &[0.3, 0.7], &[2.0], 0.0, &mut next, &mut y).unwrap();
    assert_eq!(next, [0.3, 0.7]); assert_eq!(y, [0.7]);
    assert_eq!(record.delta_h, 0.0); assert_eq!(record.supplied, 0.0); assert_eq!(record.dissipated, 0.0);
}

#[test]
fn supply_audit_still_detects_a_nonpassive_raw_system() {
    let sys = PortHamiltonian::from_raw_parts(1, 0, vec![0.2], vec![0.0], vec![],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()));
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut next = [0.0];
    let record = work.step_into(&sys, &[1.0], &[], 0.01, &mut next, &mut []).unwrap();
    assert!(record.supply_defect() > 0.001);
}

#[test]
fn cancellation_inside_jacobian_preserves_outputs_and_exact_retry() {
    let sys = nonlinear();
    let x = [0.25, 0.4, -0.13, 0.1];
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([123.0; 4], [456.0]);
    let mut polls = 0;
    let error = work.step_into_controlled(&sys, &x, &[0.2], 0.002,
        &mut next, &mut y, || { polls += 1; polls == 5 }).unwrap_err();
    assert_eq!(error, fs_phs::PreparedStepError::Cancelled);
    assert_eq!(polls, 5);
    assert_eq!(next, [123.0; 4]);
    assert_eq!(y, [456.0]);
    work.step_into(&sys, &x, &[0.2], 0.002, &mut next, &mut y).unwrap();
    let mut clean = StepWorkspace::new(&sys).unwrap();
    let (mut expected, mut expected_y) = ([0.0; 4], [0.0]);
    clean.step_into(&sys, &x, &[0.2], 0.002, &mut expected, &mut expected_y).unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(y.map(f64::to_bits), expected_y.map(f64::to_bits));
}

#[test]
fn cancellation_at_publication_does_not_leak_a_converged_candidate() {
    let sys = oscillator(0.7);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([0.0; 2], [0.0]);
    let mut total_polls = 0;
    work.step_into_controlled(&sys, &[0.2, -0.4], &[0.8], 0.003,
        &mut next, &mut y, || { total_polls += 1; false }).unwrap();
    assert!(total_polls > 5);
    next = [123.0; 2]; y = [456.0];
    let mut polls = 0;
    let error = work.step_into_controlled(&sys, &[0.2, -0.4], &[0.8], 0.003,
        &mut next, &mut y, || { polls += 1; polls == total_polls }).unwrap_err();
    assert_eq!(error, fs_phs::PreparedStepError::Cancelled);
    assert_eq!(next, [123.0; 2]);
    assert_eq!(y, [456.0]);
}

#[test]
fn prepared_workspace_reads_the_current_operators_not_a_cached_model() {
    let undamped = oscillator(0.0);
    let damped = oscillator(1.5);
    let mut work = StepWorkspace::new(&undamped).unwrap();
    let mut clean = StepWorkspace::new(&damped).unwrap();
    let (mut actual, mut actual_y) = ([0.0; 2], [0.0]);
    let (mut expected, mut expected_y) = ([0.0; 2], [0.0]);
    work.step_into(&undamped, &[0.2, 0.3], &[0.1], 0.002, &mut actual, &mut actual_y).unwrap();
    work.step_into(&damped, &[0.2, 0.3], &[0.1], 0.002, &mut actual, &mut actual_y).unwrap();
    clean.step_into(&damped, &[0.2, 0.3], &[0.1], 0.002, &mut expected, &mut expected_y).unwrap();
    assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(actual_y.map(f64::to_bits), expected_y.map(f64::to_bits));
}

#[test]
fn prepared_maxwell_history_retains_reference_storage_and_work() {
    let sys = oscillator(0.1).with_relaxation_branches(vec![fs_phs::RelaxationBranch {
        projection: vec![1.0, 0.0], stiffness: 3.0, relaxation_time_s: 0.012,
    }]).unwrap();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut x = [0.02, -0.01, 0.0];
    let mut reference = x.to_vec();
    let (mut next, mut y) = ([0.0; 3], [0.0]);
    for _ in 0..100 {
        let expected = step(&sys, &reference, &[0.02], 0.0005).unwrap();
        let record = work.step_into(&sys, &x, &[0.02], 0.0005, &mut next, &mut y).unwrap();
        assert!(record.dissipated >= 0.0);
        assert!(record.balance_residual().abs() < 1e-11);
        for (a, b) in next.iter().zip(&expected.x) { assert!((a - b).abs() < 1e-9); }
        x = next;
        reference = expected.x;
    }
    assert_ne!(x[2], 0.0);
}
