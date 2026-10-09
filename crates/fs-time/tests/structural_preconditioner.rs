//! Effective structural preconditioning, atomic failure, and checkpoint gradients.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    GeneralizedAlpha, ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderOperatorWeights,
    SecondOrderProblem, SecondOrderState, TimeSolveError, galpha_step,
    second_order_adjoint::{
        SecondOrderAdjointConfig, SecondOrderVjp,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
            StructuralReplayBudget, StructuralTrajectoryModel,
        },
    },
};
use std::cell::{Cell, RefCell};

const H: f64 = 0.05;
const RHO: f64 = 0.4;
const MASS: [f64; 4] = [1e-3, 0.0, 0.0, 1e3];
const DAMPING: [f64; 4] = [0.02, 0.003, -0.01, 2.0];
const TERMINAL: ([f64; 2], [f64; 2], [f64; 2]) = ([0.6, -0.2], [0.03, 0.02], [-0.001, 0.004]);

fn multiply(matrix: [f64; 4], input: &[f64], out: &mut [f64]) {
    out[0] = matrix[0] * input[0] + matrix[1] * input[1];
    out[1] = matrix[2] * input[0] + matrix[3] * input[1];
}

fn transpose(matrix: [f64; 4]) -> [f64; 4] {
    [matrix[0], matrix[2], matrix[1], matrix[3]]
}

fn solve(matrix: [f64; 4], rhs: &[f64], out: &mut [f64]) {
    let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    out[0] = (matrix[3] * rhs[0] - matrix[1] * rhs[1]) / determinant;
    out[1] = (matrix[0] * rhs[1] - matrix[2] * rhs[0]) / determinant;
}

fn weights() -> SecondOrderOperatorWeights {
    let am = (2.0 * RHO - 1.0) / (RHO + 1.0);
    let af = RHO / (RHO + 1.0);
    let gamma = 0.5 - am + af;
    let beta = 0.25 * (1.0 - am + af) * (1.0 - am + af);
    SecondOrderOperatorWeights {
        mass: (1.0 - am) / (beta * H * H),
        damping: (1.0 - af) * gamma / (beta * H),
        tangent: 1.0 - af,
    }
}

#[derive(Debug)]
struct Call {
    q: [f64; 2],
    weights: SecondOrderOperatorWeights,
    outer: usize,
    inner: usize,
}

struct Model {
    parameter: f64,
    // 0: exact inverse; 1: identity; 2: NaN; 3: incomplete output; 4: panic.
    policy: Cell<u8>,
    calls: RefCell<Vec<Call>>,
    tangent_points: RefCell<Vec<[f64; 2]>>,
    residual_calls: Cell<usize>,
}

impl Model {
    fn new(parameter: f64) -> Self {
        Self {
            parameter,
            policy: Cell::new(0),
            calls: RefCell::new(Vec::new()),
            tangent_points: RefCell::new(Vec::new()),
            residual_calls: Cell::new(0),
        }
    }

    fn stiffness(&self) -> [f64; 4] {
        [
            2.0 + 0.2 * self.parameter,
            0.7,
            -0.4,
            3e5 * (1.0 + 0.1 * self.parameter),
        ]
    }

    fn effective(&self, weights: SecondOrderOperatorWeights) -> [f64; 4] {
        let stiffness = self.stiffness();
        std::array::from_fn(|i| {
            weights.mass * MASS[i] + weights.damping * DAMPING[i] + weights.tangent * stiffness[i]
        })
    }
}

impl SecondOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
        multiply(MASS, x, out);
    }
    fn damping_apply(&self, x: &[f64], out: &mut [f64]) {
        multiply(DAMPING, x, out);
    }
    fn internal_force(&self, q: &[f64], out: &mut [f64]) {
        self.residual_calls.set(self.residual_calls.get() + 1);
        multiply(self.stiffness(), q, out);
    }
    fn tangent_apply(&self, q: &[f64], x: &[f64], out: &mut [f64]) {
        self.tangent_points.borrow_mut().push([q[0], q[1]]);
        multiply(self.stiffness(), x, out);
    }
    fn preconditioner_apply(
        &self,
        q: &[f64],
        weights: SecondOrderOperatorWeights,
        outer_iteration: usize,
        inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        self.calls.borrow_mut().push(Call {
            q: [q[0], q[1]],
            weights,
            outer: outer_iteration,
            inner: inner_iteration,
        });
        match self.policy.get() {
            1 => output.copy_from_slice(residual),
            2 => output.fill(f64::NAN),
            3 => output[0] = residual[0],
            4 => panic!("injected structural preconditioner failure"),
            _ => solve(self.effective(weights), residual, output),
        }
    }
}

impl SecondOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        1
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        multiply(transpose(MASS), seed, output);
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        multiply(transpose(DAMPING), seed, output);
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        _: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        multiply(transpose(self.stiffness()), seed, output);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        _: &[f64],
        _: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        output[0] = 0.2 * q[0] * seed[0] + 3e4 * q[1] * seed[1];
        Ok(())
    }
}

impl StructuralTrajectoryModel for Model {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        output[0] = 0.5 * (1.0 + 0.2 * time);
        output[1] = 2e4 * self.parameter * (1.0 + time);
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        output[0] = 2e4 * (1.0 + time) * seed[1];
        Ok(())
    }
}

struct AdjointPreconditioner([f64; 4]);
impl FlexiblePreconditioner for AdjointPreconditioner {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        solve(transpose(self.0), residual, output);
    }
}

fn initial() -> SecondOrderState {
    SecondOrderState::new(0.0, &[0.01, -0.02], &[0.02, 0.005], &[0.0, 0.0])
}

fn method() -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        2,
        H,
        RHO,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-12,
                linear_restart: 1,
                max_linear_cycles: 1,
                forcing_maximum: 1e-8,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 4,
        },
    )
}

fn dense(parameter: f64, initial: &SecondOrderState, steps: usize) -> SecondOrderState {
    let model = Model::new(parameter);
    let method = GeneralizedAlpha::new(&MASS, &DAMPING, &model.stiffness(), 2, H, RHO);
    let mut state = initial.clone();
    for _ in 0..steps {
        let time = state.t + H / (1.0 + RHO);
        let load = [0.5 * (1.0 + 0.2 * time), 2e4 * parameter * (1.0 + time)];
        galpha_step(&method, &mut state.q, &mut state.v, &mut state.a, &load);
        state.t += H;
    }
    state
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn compare(actual: &SecondOrderState, expected: &SecondOrderState) {
    for (a, b) in actual
        .q
        .iter()
        .chain(&actual.v)
        .chain(&actual.a)
        .zip(expected.q.iter().chain(&expected.v).chain(&expected.a))
    {
        close(*a, *b, 1e-10);
    }
}

#[test]
fn stage_operator_preconditioning_resolves_a_nonsymmetric_stiff_step_with_one_krylov_column() {
    let model = Model::new(0.7);
    let start = initial();
    let mut state = start.clone();
    let mut forcing = [0.0; 2];
    model
        .forcing(method().forcing_time(start.t).unwrap(), &mut forcing)
        .unwrap();
    let telemetry = method().step(&mut state, &model, &forcing).unwrap();
    assert!(telemetry.newton.converged);
    assert!(
        telemetry
            .newton
            .history
            .iter()
            .all(|iteration| iteration.linear_iterations == 1)
    );
    compare(&state, &dense(model.parameter, &start, 1));
    let calls = model.calls.borrow();
    assert!(!calls.is_empty());
    assert_eq!((calls[0].outer, calls[0].inner), (0, 0));
    for call in calls.iter() {
        assert_eq!(call.weights, weights());
        assert!(model.tangent_points.borrow().contains(&call.q));
    }
    let af = RHO / (1.0 + RHO);
    for i in 0..2 {
        let predicted = H.mul_add(start.v[i], start.q[i]);
        close(
            calls[0].q[i],
            (1.0 - af).mul_add(predicted, af * start.q[i]),
            1e-15,
        );
    }
    let identity = Model::new(0.7);
    identity.policy.set(1);
    let mut refused = start.clone();
    assert!(method().step(&mut refused, &identity, &forcing).is_err());
    assert_eq!(refused, start);
}

#[test]
fn preconditioner_faults_and_cancellation_preserve_an_existing_structural_checkpoint() {
    let model = Model::new(0.7);
    let mut state = initial();
    method().step(&mut state, &model, &[0.5, 1.4e4]).unwrap();
    let checkpoint = state.clone();
    for policy in [2, 3] {
        model.policy.set(policy);
        assert!(method().step(&mut state, &model, &[0.7, 1.5e4]).is_err());
        assert_eq!(state, checkpoint);
    }
    model.policy.set(0);
    let calls_before = model.calls.borrow().len();
    assert_eq!(
        method().step_controlled(&mut state, &model, &[0.7, 1.5e4], &mut || model
            .calls
            .borrow()
            .len()
            > calls_before),
        Err(TimeSolveError::Cancelled)
    );
    assert_eq!(state, checkpoint);
}

fn same_state_bits(actual: &SecondOrderState, expected: &SecondOrderState) {
    assert_eq!(actual, expected);
    assert_eq!(actual.t.to_bits(), expected.t.to_bits());
    assert_eq!(
        actual
            .q
            .iter()
            .chain(&actual.v)
            .chain(&actual.a)
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        expected
            .q
            .iter()
            .chain(&expected.v)
            .chain(&expected.a)
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
    );
}

#[test]
fn callback_cancellation_stops_inner_work_and_retries_the_same_structural_step() {
    let model = Model::new(0.7);
    let method = method();
    let mut checkpoint = initial();
    method.step(&mut checkpoint, &model, &[0.5, 1.4e4]).unwrap();
    let forcing = [0.7, 1.5e4];
    let mut expected = checkpoint.clone();
    method.step(&mut expected, &model, &forcing).unwrap();

    for phase in 0..3 {
        let mut state = checkpoint.clone();
        let tangent_before = model.tangent_points.borrow().len();
        let residual_before = model.residual_calls.get();
        let preconditioner_before = model.calls.borrow().len();
        let error = method
            .step_controlled(&mut state, &model, &forcing, &mut || match phase {
                0 => model.tangent_points.borrow().len() > tangent_before,
                1 => model.calls.borrow().len() > preconditioner_before,
                _ => model.residual_calls.get() > residual_before + 1,
            })
            .unwrap_err();
        assert_eq!(error, TimeSolveError::Cancelled);
        same_state_bits(&state, &checkpoint);
        let preconditioners = model.calls.borrow().len() - preconditioner_before;
        let tangents = model.tangent_points.borrow().len() - tangent_before;
        let residuals = model.residual_calls.get() - residual_before;
        match phase {
            0 => assert_eq!((preconditioners, tangents, residuals), (0, 1, 1)),
            1 => assert_eq!((preconditioners, tangents, residuals), (1, 1, 1)),
            _ => assert_eq!((preconditioners, tangents, residuals), (1, 3, 2)),
        }
        method.step(&mut state, &model, &forcing).unwrap();
        same_state_bits(&state, &expected);
    }

    let mut state = checkpoint.clone();
    model.policy.set(4);
    assert_eq!(
        method.step(&mut state, &model, &forcing),
        Err(TimeSolveError::SolverCallbackPanicked(
            fs_solver::SolverCallback::Preconditioner
        ))
    );
    same_state_bits(&state, &checkpoint);
    model.policy.set(0);
    method.step(&mut state, &model, &forcing).unwrap();
    same_state_bits(&state, &expected);
}

fn objective(state: &SecondOrderState, parameter: f64) -> f64 {
    state
        .q
        .iter()
        .zip(TERMINAL.0)
        .map(|(x, seed)| x * seed)
        .sum::<f64>()
        + state
            .v
            .iter()
            .zip(TERMINAL.1)
            .map(|(x, seed)| x * seed)
            .sum::<f64>()
        + state
            .a
            .iter()
            .zip(TERMINAL.2)
            .map(|(x, seed)| x * seed)
            .sum::<f64>()
        + 0.13 * parameter
}

#[test]
fn checkpoint_replay_uses_the_hook_and_preserves_the_dense_trajectory_gradient() {
    let model = Model::new(0.7);
    let start = initial();
    let mut recorded = RecordedStructural::new(
        method(),
        &model,
        &start,
        StructuralRecordingConfig {
            steps: 5,
            adjoint: SecondOrderAdjointConfig {
                restart: 1,
                max_cycles: 2,
                tolerance: 1e-12,
            },
            max_workspace_components: 2000,
        },
    )
    .unwrap();
    assert_eq!(
        recorded.advance(5, 5, &mut || false).unwrap().status,
        StructuralRecordingStatus::ReachedEnd
    );
    compare(recorded.state(), &dense(model.parameter, &start, 5));
    let calls_before = model.calls.borrow().len();
    let gradient = recorded
        .pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &[0.13],
            &AdjointPreconditioner(model.effective(weights())),
            StructuralReplayBudget {
                checkpoints: 4,
                forward_steps: 32,
            },
            &mut || false,
        )
        .unwrap();
    assert!(model.calls.borrow().len() > calls_before);
    assert!(gradient.replayed_steps >= 5);
    let epsilon = 1e-5;
    let fd = (objective(
        &dense(model.parameter + epsilon, &start, 5),
        model.parameter + epsilon,
    ) - objective(
        &dense(model.parameter - epsilon, &start, 5),
        model.parameter - epsilon,
    )) / (2.0 * epsilon);
    close(gradient.parameters[0], fd, 2e-7);
    for (component, derivative) in [
        &gradient.initial_q,
        &gradient.initial_v,
        &gradient.initial_a,
    ]
    .into_iter()
    .enumerate()
    {
        for i in 0..2 {
            let mut plus = start.clone();
            let mut minus = start.clone();
            [&mut plus.q, &mut plus.v, &mut plus.a][component][i] += epsilon;
            [&mut minus.q, &mut minus.v, &mut minus.a][component][i] -= epsilon;
            let fd = (objective(&dense(model.parameter, &plus, 5), model.parameter)
                - objective(&dense(model.parameter, &minus, 5), model.parameter))
                / (2.0 * epsilon);
            close(derivative[i], fd, 2e-7);
        }
    }
}
