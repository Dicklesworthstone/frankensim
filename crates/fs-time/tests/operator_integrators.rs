//! G1/G3 evidence for operator-backed implicit time integrators.

use fs_solver::{Globalization, LineSearchConfig, LinearOp, NewtonKrylovConfig};
use fs_time::{
    FirstOrderGeneralizedAlpha, FirstOrderProblem, FirstOrderState, GeneralizedAlpha,
    IdentityPreconditioner, Imex2, ImexSolveConfig, ImexState, ImplicitSolveConfig,
    LinearFirstOrderSystem, LinearSecondOrderSystem, OperatorFirstOrderGeneralizedAlpha,
    OperatorGeneralizedAlpha, OperatorImex2, SecondOrderState, first_order_galpha_step,
    galpha_step, imex2_step,
};

#[derive(Debug, Clone)]
struct DenseOp {
    n: usize,
    entries: Vec<f64>,
}

impl DenseOp {
    fn new(n: usize, entries: &[f64]) -> Self {
        assert_eq!(entries.len(), n * n);
        Self {
            n,
            entries: entries.to_vec(),
        }
    }
}

impl LinearOp for DenseOp {
    fn n(&self) -> usize {
        self.n
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) {
        for (row, value) in output.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (column, input_value) in input.iter().copied().enumerate() {
                sum = self.entries[row * self.n + column].mul_add(input_value, sum);
            }
            *value = sum;
        }
    }

    fn apply_transpose(&self, input: &[f64], output: &mut [f64]) {
        for (column, value) in output.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (row, input_value) in input.iter().copied().enumerate() {
                sum = self.entries[row * self.n + column].mul_add(input_value, sum);
            }
            *value = sum;
        }
    }
}

#[derive(Debug, Clone)]
struct DiagonalOp(Vec<f64>);

impl LinearOp for DiagonalOp {
    fn n(&self) -> usize {
        self.0.len()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) {
        for ((output, diagonal), input) in output.iter_mut().zip(&self.0).zip(input) {
            *output = *diagonal * *input;
        }
    }
}

fn tight_newton() -> ImplicitSolveConfig {
    ImplicitSolveConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1.0e-11,
            relative_tolerance: 1.0e-11,
            linear_restart: 8,
            max_linear_cycles: 4,
            forcing_minimum: 1.0e-12,
            forcing_maximum: 1.0e-12,
            forcing_gamma: 0.9,
            forcing_exponent: 1.5,
            globalization: Globalization::LineSearch(LineSearchConfig::default()),
        },
        max_newton_iterations: 8,
    }
}

fn imex_config() -> ImexSolveConfig {
    ImexSolveConfig {
        tolerance: 1.0e-12,
        restart: 8,
        max_cycles: 4,
    }
}

fn assert_close(left: &[f64], right: &[f64], tolerance: f64) {
    assert_eq!(left.len(), right.len());
    for (index, (left, right)) in left.iter().zip(right).enumerate() {
        assert!(
            (left - right).abs() <= tolerance,
            "entry {index}: {left:.16e} != {right:.16e} (tol {tolerance:.3e})"
        );
    }
}

fn assert_same_bits(left: &[f64], right: &[f64]) {
    assert_eq!(left.len(), right.len());
    for (index, (left, right)) in left.iter().zip(right).enumerate() {
        assert_eq!(
            left.to_bits(),
            right.to_bits(),
            "entry {index} differs bitwise"
        );
    }
}

#[test]
fn structural_dense_operator_agreement_and_split_replay() {
    let mass = [1.0, 0.0, 0.0, 1.0];
    let damping = [0.04, 0.01, 0.01, 0.06];
    let stiffness = [4.0, -0.2, -0.2, 9.0];
    let h = 0.025;
    let rho = 0.7;
    let dense = GeneralizedAlpha::new(&mass, &damping, &stiffness, 2, h, rho);
    let mass_op = DenseOp::new(2, &mass);
    let damping_op = DenseOp::new(2, &damping);
    let stiffness_op = DenseOp::new(2, &stiffness);
    let system = LinearSecondOrderSystem::new(&mass_op, &damping_op, &stiffness_op);
    let operator = OperatorGeneralizedAlpha::new(2, h, rho, tight_newton());
    let q0 = [1.0, 0.5];
    let v0 = [0.0, 0.1];
    let a0 = [-3.9, -4.306];
    let forcing = [0.0, 0.0];

    let (mut q_dense, mut v_dense, mut a_dense) = (q0.to_vec(), v0.to_vec(), a0.to_vec());
    let mut straight = SecondOrderState::new(0.0, &q0, &v0, &a0);
    for _ in 0..40 {
        galpha_step(&dense, &mut q_dense, &mut v_dense, &mut a_dense, &forcing);
        operator.step(&mut straight, &system, &forcing).unwrap();
    }
    assert_close(&q_dense, &straight.q, 2.0e-10);
    assert_close(&v_dense, &straight.v, 2.0e-10);
    assert_close(&a_dense, &straight.a, 2.0e-9);

    let mut prefix = SecondOrderState::new(0.0, &q0, &v0, &a0);
    for _ in 0..13 {
        operator.step(&mut prefix, &system, &forcing).unwrap();
    }
    let mut resumed = prefix.clone();
    for _ in 13..40 {
        operator.step(&mut resumed, &system, &forcing).unwrap();
    }
    assert_same_bits(&straight.q, &resumed.q);
    assert_same_bits(&straight.v, &resumed.v);
    assert_same_bits(&straight.a, &resumed.a);
    assert_eq!(straight.t.to_bits(), resumed.t.to_bits());
    assert_eq!(straight.steps, resumed.steps);
    assert_eq!(straight.history, resumed.history);
    let krylov: usize = straight
        .history
        .iter()
        .map(|step| step.krylov_iterations())
        .sum();
    println!(
        "{{\"case\":\"structural-dense-operator\",\"steps\":{},\"newton\":{},\"krylov\":{},\"q0\":{:.12e}}}",
        straight.steps,
        straight
            .history
            .iter()
            .map(|step| step.newton.iterations)
            .sum::<usize>(),
        krylov,
        straight.q[0]
    );
}

fn structural_error(h: f64, operator_path: bool) -> f64 {
    let steps = (1.0 / h).round() as usize;
    let mass = [1.0];
    let damping = [0.0];
    let stiffness = [1.0];
    if operator_path {
        let mass_op = DiagonalOp(mass.to_vec());
        let damping_op = DiagonalOp(damping.to_vec());
        let stiffness_op = DiagonalOp(stiffness.to_vec());
        let system = LinearSecondOrderSystem::new(&mass_op, &damping_op, &stiffness_op);
        let method = OperatorGeneralizedAlpha::new(1, h, 0.8, tight_newton());
        let mut state = SecondOrderState::new(0.0, &[1.0], &[0.0], &[-1.0]);
        for _ in 0..steps {
            method.step(&mut state, &system, &[0.0]).unwrap();
        }
        (state.q[0] - 1.0f64.cos())
            .abs()
            .max((state.v[0] + 1.0f64.sin()).abs())
    } else {
        let method = GeneralizedAlpha::new(&mass, &damping, &stiffness, 1, h, 0.8);
        let (mut q, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
        for _ in 0..steps {
            galpha_step(&method, &mut q, &mut v, &mut a, &[0.0]);
        }
        (q[0] - 1.0f64.cos()).abs().max((v[0] + 1.0f64.sin()).abs())
    }
}

#[test]
fn structural_generalized_alpha_is_second_order_on_both_paths() {
    for operator_path in [false, true] {
        let coarse = structural_error(0.1, operator_path);
        let fine = structural_error(0.05, operator_path);
        let order = (coarse / fine).log2();
        assert!(
            (order - 2.0).abs() < 0.35,
            "operator={operator_path} order={order:.3} errors={coarse:.3e}/{fine:.3e}"
        );
        println!(
            "{{\"case\":\"structural-order\",\"operator\":{operator_path},\"order\":{order:.6},\"coarse\":{coarse:.12e},\"fine\":{fine:.12e}}}"
        );
    }
}

fn first_order_linear_error(h: f64, operator_path: bool) -> f64 {
    let steps = (1.0 / h).round() as usize;
    if operator_path {
        let mass = DiagonalOp(vec![1.0]);
        let evolution = DiagonalOp(vec![1.0]);
        let system = LinearFirstOrderSystem::new(&mass, &evolution);
        let method = OperatorFirstOrderGeneralizedAlpha::new(1, h, 0.6, tight_newton());
        let mut state = FirstOrderState::new(0.0, &[1.0], &[-1.0]);
        for _ in 0..steps {
            method.step(&mut state, &system, &[0.0]).unwrap();
        }
        (state.u[0] - (-1.0f64).exp()).abs()
    } else {
        let method = FirstOrderGeneralizedAlpha::new(&[1.0], &[1.0], 1, h, 0.6);
        let (mut u, mut rate) = (vec![1.0], vec![-1.0]);
        for _ in 0..steps {
            first_order_galpha_step(&method, &mut u, &mut rate, &[0.0]);
        }
        (u[0] - (-1.0f64).exp()).abs()
    }
}

#[test]
fn first_order_generalized_alpha_dense_operator_order_and_agreement() {
    for operator_path in [false, true] {
        let coarse = first_order_linear_error(0.1, operator_path);
        let fine = first_order_linear_error(0.05, operator_path);
        let order = (coarse / fine).log2();
        assert!(
            (order - 2.0).abs() < 0.35,
            "operator={operator_path} first-order-system order={order:.3}"
        );
    }

    let h = 0.025;
    let dense = FirstOrderGeneralizedAlpha::new(&[2.0], &[3.0], 1, h, 0.4);
    let mass = DenseOp::new(1, &[2.0]);
    let evolution = DenseOp::new(1, &[3.0]);
    let system = LinearFirstOrderSystem::new(&mass, &evolution);
    let operator = OperatorFirstOrderGeneralizedAlpha::new(1, h, 0.4, tight_newton());
    let (mut u_dense, mut rate_dense) = (vec![0.75], vec![-1.125]);
    let mut state = FirstOrderState::new(0.0, &u_dense, &rate_dense);
    for _ in 0..40 {
        first_order_galpha_step(&dense, &mut u_dense, &mut rate_dense, &[0.0]);
        operator.step(&mut state, &system, &[0.0]).unwrap();
    }
    assert_close(&u_dense, &state.u, 2.0e-11);
    assert_close(&rate_dense, &state.rate, 2.0e-10);

    let mut prefix = FirstOrderState::new(0.0, &[0.75], &[-1.125]);
    for _ in 0..15 {
        operator.step(&mut prefix, &system, &[0.0]).unwrap();
    }
    let mut resumed = prefix.clone();
    for _ in 15..40 {
        operator.step(&mut resumed, &system, &[0.0]).unwrap();
    }
    assert_same_bits(&state.u, &resumed.u);
    assert_same_bits(&state.rate, &resumed.rate);
    assert_eq!(state.t.to_bits(), resumed.t.to_bits());
    assert_eq!(state.history, resumed.history);
    println!(
        "{{\"case\":\"first-order-dense-operator\",\"steps\":{},\"u\":{:.12e},\"rate\":{:.12e}}}",
        state.steps, state.u[0], state.rate[0]
    );
}

struct CubicDecay;

impl FirstOrderProblem for CubicDecay {
    fn dimension(&self) -> usize {
        1
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        output[0] = input[0];
    }

    fn internal_force(&self, _t: f64, u: &[f64], output: &mut [f64]) {
        output[0] = 0.1f64.mul_add(u[0] * u[0] * u[0], u[0]);
    }

    fn tangent_apply(&self, _t: f64, u: &[f64], direction: &[f64], output: &mut [f64]) {
        output[0] = (1.0 + 0.3 * u[0] * u[0]) * direction[0];
    }
}

fn nonlinear_first_order_endpoint(h: f64) -> f64 {
    let steps = (1.0 / h).round() as usize;
    let method = OperatorFirstOrderGeneralizedAlpha::new(1, h, 1.0, tight_newton());
    let mut state = FirstOrderState::new(0.0, &[0.8], &[-0.8512]);
    for _ in 0..steps {
        method.step(&mut state, &CubicDecay, &[0.0]).unwrap();
    }
    state.u[0]
}

#[test]
fn nonlinear_first_order_self_convergence_is_second_order() {
    let (coarse, medium, fine) = (
        nonlinear_first_order_endpoint(0.1),
        nonlinear_first_order_endpoint(0.05),
        nonlinear_first_order_endpoint(0.025),
    );
    let coarse_delta = (coarse - medium).abs();
    let fine_delta = (medium - fine).abs();
    let order = (coarse_delta / fine_delta).log2();
    assert!(
        (order - 2.0).abs() < 0.4,
        "nonlinear self-convergence order={order:.3}, deltas={coarse_delta:.3e}/{fine_delta:.3e}"
    );
    println!(
        "{{\"case\":\"first-order-nonlinear-order\",\"order\":{order:.6},\"coarse_delta\":{coarse_delta:.12e},\"fine_delta\":{fine_delta:.12e}}}"
    );
}

#[test]
fn imex_dense_operator_agreement_order_and_split_replay() {
    let h = 0.02;
    let dense = Imex2::new(&[-2.0, 0.25, 0.0, -0.5], 2, h);
    let linear = DenseOp::new(2, &[-2.0, 0.25, 0.0, -0.5]);
    let operator = OperatorImex2::new(2, h, imex_config());
    let nonlin = |u: &[f64], output: &mut [f64]| {
        output[0] = u[0] * u[0];
        output[1] = -0.2 * u[1] * u[1] * u[1];
    };
    let mut u_dense = vec![0.5, -0.25];
    let mut straight = ImexState::new(0.0, &u_dense);
    for _ in 0..50 {
        imex2_step(&dense, &mut u_dense, &nonlin);
        operator
            .step(&mut straight, &linear, &IdentityPreconditioner, &nonlin)
            .unwrap();
    }
    assert_close(&u_dense, &straight.u, 2.0e-10);

    let mut prefix = ImexState::new(0.0, &[0.5, -0.25]);
    for _ in 0..17 {
        operator
            .step(&mut prefix, &linear, &IdentityPreconditioner, &nonlin)
            .unwrap();
    }
    let mut resumed = prefix.clone();
    for _ in 17..50 {
        operator
            .step(&mut resumed, &linear, &IdentityPreconditioner, &nonlin)
            .unwrap();
    }
    assert_same_bits(&straight.u, &resumed.u);
    assert_eq!(straight.t.to_bits(), resumed.t.to_bits());
    assert_eq!(straight.steps, resumed.steps);
    assert_eq!(straight.history, resumed.history);

    let error = |step: f64| -> f64 {
        let method = OperatorImex2::new(1, step, imex_config());
        let linear = DiagonalOp(vec![-1.0]);
        let nonlinear = |u: &[f64], output: &mut [f64]| output[0] = u[0] * u[0];
        let mut state = ImexState::new(0.0, &[0.5]);
        for _ in 0..(1.0 / step).round() as usize {
            method
                .step(&mut state, &linear, &IdentityPreconditioner, &nonlinear)
                .unwrap();
        }
        let exact = 0.5 / (0.5 + 0.5 * 1.0f64.exp());
        (state.u[0] - exact).abs()
    };
    let (coarse, fine) = (error(0.1), error(0.05));
    let order = (coarse / fine).log2();
    assert!((order - 2.0).abs() < 0.35, "operator IMEX order={order:.3}");
    println!(
        "{{\"case\":\"imex-dense-operator\",\"steps\":{},\"order\":{order:.6},\"stage1_iters\":{},\"stage2_iters\":{}}}",
        straight.steps,
        straight
            .history
            .iter()
            .map(|row| row.stage_one.iters)
            .sum::<usize>(),
        straight
            .history
            .iter()
            .map(|row| row.stage_two.iters)
            .sum::<usize>()
    );
}
