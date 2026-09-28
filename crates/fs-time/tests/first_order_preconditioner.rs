//! G3/G4: the first-order Newton preconditioner uses the exact stage operator.
use fs_solver::{NewtonKrylovConfig, NewtonStallDiagnosis, StallDiagnosis};
use fs_time::galpha::{
    FirstOrderGeneralizedAlpha, FirstOrderOperatorWeights, FirstOrderProblem, FirstOrderState,
    ImplicitSolveConfig, OperatorFirstOrderGeneralizedAlpha, TimeSolveError,
    first_order_galpha_step,
};
use std::cell::{Cell, RefCell};

const MASS: [f64; 4] = [0.001, 0.0002, -0.03, 10.0];
const OPERATOR: [f64; 4] = [1e6, 25.0, -7.0, 2.0];
const H: f64 = 0.08;
const FORCING: [f64; 2] = [1.0, -2.0];

#[derive(Clone, Copy)]
enum Policy {
    Exact,
    Identity,
    NonFinite,
    Unwritten,
}

#[derive(Clone, Debug)]
struct Call {
    time: f64,
    state: [f64; 2],
    weights: [f64; 2],
    outer: usize,
    inner: usize,
    tangent: (f64, [f64; 2]),
}

struct Model {
    policy: Cell<Policy>,
    nonlinear: bool,
    calls: RefCell<Vec<Call>>,
    tangent: Cell<(f64, [f64; 2])>,
}

impl Model {
    fn new(policy: Policy, nonlinear: bool) -> Self {
        Self {
            policy: Cell::new(policy),
            nonlinear,
            calls: RefCell::new(Vec::new()),
            tangent: Cell::new((f64::NAN, [f64::NAN; 2])),
        }
    }
    fn jacobian(&self, u: &[f64]) -> [f64; 4] {
        let mut matrix = OPERATOR;
        if self.nonlinear {
            matrix[3] += 360.0 * u[1] * u[1];
        }
        matrix
    }
}

fn apply(matrix: [f64; 4], input: &[f64], out: &mut [f64]) {
    out[0] = matrix[0] * input[0] + matrix[1] * input[1];
    out[1] = matrix[2] * input[0] + matrix[3] * input[1];
}

impl FirstOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, input: &[f64], out: &mut [f64]) {
        apply(MASS, input, out);
    }
    fn internal_force(&self, _: f64, u: &[f64], out: &mut [f64]) {
        apply(OPERATOR, u, out);
        if self.nonlinear {
            out[1] += 120.0 * u[1].powi(3);
        }
    }
    fn tangent_apply(&self, time: f64, u: &[f64], input: &[f64], out: &mut [f64]) {
        self.tangent.set((time, [u[0], u[1]]));
        apply(self.jacobian(u), input, out);
    }
    fn preconditioner_apply(
        &self,
        time: f64,
        u: &[f64],
        weights: FirstOrderOperatorWeights,
        outer_iteration: usize,
        inner_iteration: usize,
        residual: &[f64],
        out: &mut [f64],
    ) {
        self.calls.borrow_mut().push(Call {
            time,
            state: [u[0], u[1]],
            weights: [weights.mass, weights.tangent],
            outer: outer_iteration,
            inner: inner_iteration,
            tangent: self.tangent.get(),
        });
        match self.policy.get() {
            Policy::Identity => out.copy_from_slice(residual),
            Policy::NonFinite => out.fill(f64::NAN),
            Policy::Unwritten => out[0] = residual[0],
            Policy::Exact => {
                let tangent = self.jacobian(u);
                let a: [f64; 4] =
                    std::array::from_fn(|i| weights.mass * MASS[i] + weights.tangent * tangent[i]);
                let determinant = a[0] * a[3] - a[1] * a[2];
                out[0] = (a[3] * residual[0] - a[1] * residual[1]) / determinant;
                out[1] = (a[0] * residual[1] - a[2] * residual[0]) / determinant;
            }
        }
    }
}

/// Implements only the original model contract so identity is truly omitted.
struct DefaultIdentity<'a>(&'a Model);
impl FirstOrderProblem for DefaultIdentity<'_> {
    fn dimension(&self) -> usize {
        self.0.dimension()
    }
    fn mass_apply(&self, input: &[f64], out: &mut [f64]) {
        self.0.mass_apply(input, out);
    }
    fn internal_force(&self, time: f64, u: &[f64], out: &mut [f64]) {
        self.0.internal_force(time, u, out);
    }
    fn tangent_apply(&self, time: f64, u: &[f64], input: &[f64], out: &mut [f64]) {
        self.0.tangent_apply(time, u, input, out);
    }
}

fn method(rho: f64, restart: usize, cycles: usize) -> OperatorFirstOrderGeneralizedAlpha {
    OperatorFirstOrderGeneralizedAlpha::new(
        2,
        H,
        rho,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-13,
                linear_restart: restart,
                max_linear_cycles: cycles,
                forcing_minimum: 1e-14,
                forcing_maximum: 1e-9,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 12,
        },
    )
}

fn initial() -> FirstOrderState {
    FirstOrderState::new(0.7, &[1.2e-5, -0.3], &[0.001, 0.1])
}

fn same_state(actual: &FirstOrderState, expected: &FirstOrderState) {
    assert_eq!(actual, expected);
    assert_eq!(actual.t.to_bits(), expected.t.to_bits());
    assert_eq!(
        actual
            .u
            .iter()
            .chain(&actual.rate)
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        expected
            .u
            .iter()
            .chain(&expected.rate)
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
    );
}

#[test]
fn a_stiff_nonsymmetric_step_converges_with_one_preconditioned_krylov_column() {
    let model = Model::new(Policy::Exact, false);
    let method = method(0.35, 1, 1);
    let mut state = initial();
    let mut refused = state.clone();
    assert!(matches!(
        method.step(&mut refused, &DefaultIdentity(&model), &FORCING),
        Err(TimeSolveError::NotConverged(_))
    ));
    same_state(&refused, &state);
    let dense = FirstOrderGeneralizedAlpha::new(&MASS, &OPERATOR, 2, H, 0.35);
    let mut u = state.u.clone();
    let mut rate = state.rate.clone();
    first_order_galpha_step(&dense, &mut u, &mut rate, &FORCING);
    let report = method.step(&mut state, &model, &FORCING).unwrap();
    assert!(report.newton.converged);
    assert!(
        report
            .newton
            .history
            .iter()
            .all(|iteration| iteration.linear_iterations == 1)
    );
    for (actual, expected) in state.u.iter().chain(&state.rate).zip(u.iter().chain(&rate)) {
        assert!(
            (actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
            "{actual} != {expected}"
        );
    }
    assert_eq!(state.steps, 1);
    assert_eq!(state.history.len(), 1);
}

#[test]
fn preconditioner_receives_the_actual_nonlinear_stage_and_logical_iteration_keys() {
    let rho = 0.35;
    let model = Model::new(Policy::Exact, true);
    let original = initial();
    let mut state = original.clone();
    method(rho, 1, 1)
        .step(&mut state, &model, &FORCING)
        .unwrap();
    let alpha_m = (3.0 - rho) / (2.0 * (1.0 + rho));
    let alpha_f = 1.0 / (1.0 + rho);
    let gamma = 0.5 + alpha_m - alpha_f;
    let calls = model.calls.borrow();
    assert!(calls.len() > 1);
    for (index, call) in calls.iter().enumerate() {
        assert_eq!(
            call.time.to_bits(),
            H.mul_add(alpha_f, original.t).to_bits()
        );
        assert_eq!(call.weights[0], alpha_m / (gamma * H));
        assert_eq!(call.weights[1], alpha_f);
        assert_eq!((call.time, call.state), call.tangent);
        assert_eq!(call.outer, index);
        assert_eq!(call.inner, 0);
    }
    for i in 0..2 {
        let guess = H.mul_add(original.rate[i], original.u[i]);
        let expected = alpha_f.mul_add(guess, (1.0 - alpha_f) * original.u[i]);
        assert_eq!(calls[0].state[i].to_bits(), expected.to_bits());
    }
    assert_ne!(calls[0].state, calls[1].state);
}

#[test]
fn faulty_preconditioners_and_cancellation_preserve_a_complete_prior_step() {
    let model = Model::new(Policy::Exact, true);
    let method = method(0.35, 1, 1);
    let mut state = initial();
    method.step(&mut state, &model, &FORCING).unwrap();
    let before = state.clone();
    assert_eq!(before.history.len(), 1);
    for policy in [Policy::NonFinite, Policy::Unwritten] {
        model.policy.set(policy);
        let error = method.step(&mut state, &model, &FORCING).unwrap_err();
        let TimeSolveError::NotConverged(report) = error else {
            panic!("unexpected error: {error:?}");
        };
        assert_eq!(
            report.diagnosis,
            Some(NewtonStallDiagnosis::LinearSolveFailed(
                StallDiagnosis::Breakdown
            ))
        );
        same_state(&state, &before);
    }
    model.policy.set(Policy::Exact);
    let calls = model.calls.borrow().len();
    assert!(matches!(
        method.step_controlled(&mut state, &model, &FORCING, &mut || model
            .calls
            .borrow()
            .len()
            > calls),
        Err(TimeSolveError::Cancelled)
    ));
    same_state(&state, &before);
    method.step(&mut state, &model, &FORCING).unwrap();
    assert_eq!(state.steps, before.steps + 1);
}

#[test]
fn explicit_identity_and_the_omitted_default_preserve_state_and_history_bits() {
    for rho in [0.0, 0.35, 1.0] {
        let model = Model::new(Policy::Identity, true);
        let mut explicit = initial();
        let mut implicit = initial();
        let method = method(rho, 2, 4);
        for _ in 0..4 {
            method.step(&mut explicit, &model, &FORCING).unwrap();
            method
                .step(&mut implicit, &DefaultIdentity(&model), &FORCING)
                .unwrap();
            same_state(&explicit, &implicit);
        }
        assert!(model.calls.borrow().iter().any(|call| call.inner > 0));
    }
}
