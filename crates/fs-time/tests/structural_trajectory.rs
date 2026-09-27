//! G3/G4/G5: structural checkpoint replay, load derivatives, and bounded resume.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    TimeSolveError,
    second_order_adjoint::{
        SecondOrderAdjointConfig, SecondOrderAdjointError, SecondOrderVjp,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
            StructuralReplayBudget, StructuralTrajectoryError, StructuralTrajectoryModel,
        },
    },
};
use std::cell::Cell;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Fault {
    None,
    ChangedLoad,
    LoadError,
    LoadIncomplete,
    DerivativeError,
    DerivativeIncomplete,
}

struct Model {
    p: [f64; 2],
    fault: Cell<Fault>,
    residual_calls: Cell<usize>,
    load_calls: Cell<usize>,
    load_vjp_calls: Cell<usize>,
}

impl Model {
    fn new(p: [f64; 2]) -> Self {
        Self {
            p,
            fault: Cell::new(Fault::None),
            residual_calls: Cell::new(0),
            load_calls: Cell::new(0),
            load_vjp_calls: Cell::new(0),
        }
    }

    fn tangent(&self, q: &[f64]) -> [f64; 4] {
        [
            2.0 + self.p[0] + 0.24 * q[0] * q[0],
            0.35,
            -0.2,
            1.5 + 0.2 * self.p[0] + 0.12 * q[1] * q[1],
        ]
    }

    fn damping(&self) -> [f64; 4] {
        [0.25 + 0.1 * self.p[1], 0.04, -0.03, 0.2 + 0.05 * self.p[1]]
    }
}

const MASS: [f64; 4] = [1.5, 0.16, -0.04, 1.2];

fn apply(a: [f64; 4], x: &[f64], y: &mut [f64]) {
    y[0] = a[0] * x[0] + a[1] * x[1];
    y[1] = a[2] * x[0] + a[3] * x[1];
}

fn transpose(a: [f64; 4], x: &[f64], y: &mut [f64]) {
    apply([a[0], a[2], a[1], a[3]], x, y);
}

impl SecondOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        apply(MASS, input, output);
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        apply(self.damping(), input, output);
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        self.residual_calls.set(self.residual_calls.get() + 1);
        output[0] = (2.0 + self.p[0]) * q[0] + 0.35 * q[1] + 0.08 * q[0].powi(3);
        output[1] = -0.2 * q[0] + (1.5 + 0.2 * self.p[0]) * q[1] + 0.04 * q[1].powi(3);
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        apply(self.tangent(q), direction, output);
    }
}

impl SecondOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        transpose(MASS, seed, output);
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        transpose(self.damping(), seed, output);
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        transpose(self.tangent(q), seed, output);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        _: &[f64],
        seed: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        out[0] = q[0] * seed[0] + 0.2 * q[1] * seed[1];
        out[1] = 0.1 * v[0] * seed[0] + 0.05 * v[1] * seed[1];
        Ok(())
    }
}

impl StructuralTrajectoryModel for Model {
    fn forcing(&self, time: f64, out: &mut [f64]) -> Result<(), String> {
        self.load_calls.set(self.load_calls.get() + 1);
        if self.fault.get() == Fault::LoadError {
            return Err("load unavailable".into());
        }
        out[0] = self.p[1] * (1.0 + 0.2 * time);
        if self.fault.get() == Fault::ChangedLoad {
            out[0] += 0.1;
        }
        if self.fault.get() == Fault::LoadIncomplete {
            return Ok(());
        }
        out[1] = 0.3 + 0.1 * time * time * self.p[1] + 0.05 * time * self.p[0];
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.load_vjp_calls.set(self.load_vjp_calls.get() + 1);
        if self.fault.get() == Fault::DerivativeError {
            return Err("load derivative unavailable".into());
        }
        out[0] = 0.05 * time * seed[1];
        if self.fault.get() == Fault::DerivativeIncomplete {
            return Ok(());
        }
        out[1] = (1.0 + 0.2 * time) * seed[0] + 0.1 * time * time * seed[1];
        Ok(())
    }
}

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, residual: &[f64], out: &mut [f64]) {
        out.copy_from_slice(residual);
    }
}

const H: f64 = 0.15;
const RHO: f64 = 0.35;
const TERMINAL: ([f64; 2], [f64; 2], [f64; 2]) = ([0.7, -0.2], [0.1, 0.3], [-0.04, 0.08]);
const DIRECT: [f64; 2] = [0.13, -0.07];

fn method() -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        2,
        H,
        RHO,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-12,
                relative_tolerance: 2e-13,
                linear_restart: 2,
                max_linear_cycles: 4,
                forcing_maximum: 1e-3,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 12,
        },
    )
}

fn config(steps: usize) -> StructuralRecordingConfig {
    let adjoint = SecondOrderAdjointConfig {
        restart: 2,
        max_cycles: 4,
        tolerance: 1e-13,
    };
    StructuralRecordingConfig {
        steps,
        adjoint,
        max_workspace_components: method().adjoint_workspace_components(2, adjoint).unwrap(),
    }
}

fn initial() -> SecondOrderState {
    let mut state = SecondOrderState::new(0.3, &[0.8, -0.3], &[0.2, 0.1], &[-0.1, 0.2]);
    state.steps = 11;
    state
}

fn recording(model: &Model, steps: usize) -> RecordedStructural<'_, Model> {
    RecordedStructural::new(method(), model, &initial(), config(steps)).unwrap()
}

fn same_state(a: &SecondOrderState, b: &SecondOrderState) {
    assert_eq!(a.t.to_bits(), b.t.to_bits());
    assert_eq!(a.steps, b.steps);
    assert_eq!(
        a.q.iter()
            .chain(&a.v)
            .chain(&a.a)
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        b.q.iter()
            .chain(&b.v)
            .chain(&b.a)
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
    );
}

fn load(model: &Model, time: f64) -> [f64; 2] {
    let stage = H.mul_add(1.0 - RHO / (1.0 + RHO), time);
    assert_eq!(
        method().forcing_time(time).unwrap().to_bits(),
        stage.to_bits()
    );
    let mut forcing = [f64::NAN; 2];
    model.forcing(stage, &mut forcing).unwrap();
    forcing
}

fn full_storage(model: &Model, steps: usize) -> (Vec<SecondOrderState>, Vec<usize>) {
    let mut states = vec![initial()];
    let mut costs = Vec::new();
    for _ in 0..steps {
        let mut state = states.last().unwrap().clone();
        let forcing = load(model, state.t);
        let before = model.residual_calls.get();
        method().step(&mut state, model, &forcing).unwrap();
        costs.push(model.residual_calls.get() - before);
        states.push(state);
    }
    (states, costs)
}

// Independent accounting of the binary schedule: replay the left prefix to
// establish the midpoint, reverse the right half, then reverse the left half.
fn replay_costs(costs: &[usize]) -> (usize, usize) {
    if costs.len() <= 1 {
        return (costs.len(), costs.iter().sum());
    }
    let mid = costs.len().div_ceil(2);
    let left = replay_costs(&costs[..mid]);
    let right = replay_costs(&costs[mid..]);
    (
        mid + left.0 + right.0,
        costs[..mid].iter().sum::<usize>() + left.1 + right.1,
    )
}

fn budget(steps: usize) -> StructuralReplayBudget {
    StructuralReplayBudget {
        checkpoints: steps.bit_width() as usize,
        forward_steps: replay_costs(&vec![0; steps]).0,
    }
}

fn objective(p: [f64; 2], steps: usize) -> f64 {
    let (states, _) = full_storage(&Model::new(p), steps);
    let end = states.last().unwrap();
    end.q
        .iter()
        .chain(&end.v)
        .chain(&end.a)
        .zip(TERMINAL.0.iter().chain(&TERMINAL.1).chain(&TERMINAL.2))
        .map(|(value, seed)| value * seed)
        .sum::<f64>()
        + p.iter()
            .zip(DIRECT)
            .map(|(value, seed)| value * seed)
            .sum::<f64>()
}

fn assert_fault(error: &StructuralTrajectoryError, fault: Fault) {
    let expected = match fault {
        Fault::LoadError => matches!(error, StructuralTrajectoryError::Forcing(_)),
        Fault::LoadIncomplete => matches!(error, StructuralTrajectoryError::NonFiniteForcing),
        Fault::DerivativeError => matches!(error, StructuralTrajectoryError::ForcingDerivative(_)),
        Fault::DerivativeIncomplete => matches!(
            error,
            StructuralTrajectoryError::Step(SecondOrderAdjointError::NonFiniteDerivative)
        ),
        _ => false,
    };
    assert!(expected, "unexpected refusal: {error:?}");
}

#[test]
fn checkpointed_loaded_trajectories_match_full_storage_and_parameter_differences() {
    for steps in [0, 1, 7, 64] {
        let model = Model::new([1.3, 0.6]);
        let mut recorded = recording(&model, steps);
        let report = recorded.advance(steps, steps, &mut || false).unwrap();
        assert_eq!(report.status, StructuralRecordingStatus::ReachedEnd);
        assert_eq!(report.advanced, steps);
        assert_eq!(recorded.accepted_steps(), steps);
        assert_eq!(recorded.state().history, Vec::new());
        let (states, costs) = full_storage(&model, steps);
        same_state(recorded.state(), states.last().unwrap());
        assert_eq!(
            recorded.time().to_bits(),
            states.last().unwrap().t.to_bits()
        );
        let (mut q, mut v, mut a) = (
            TERMINAL.0.to_vec(),
            TERMINAL.1.to_vec(),
            TERMINAL.2.to_vec(),
        );
        let mut parameters = DIRECT.to_vec();
        for state in states[..steps].iter().rev() {
            let gradient = method()
                .step_vjp(
                    state,
                    &model,
                    &load(&model, state.t),
                    (&q, &v, &a),
                    &Identity,
                    config(steps).adjoint,
                    config(steps).max_workspace_components,
                    &mut || false,
                )
                .unwrap();
            let mut load_parameters = [f64::NAN; 2];
            model
                .forcing_vjp(
                    method().forcing_time(state.t).unwrap(),
                    &gradient.forcing,
                    &mut load_parameters,
                )
                .unwrap();
            for i in 0..2 {
                parameters[i] += gradient.parameters[i];
                parameters[i] += load_parameters[i];
            }
            q = gradient.initial_q;
            v = gradient.initial_v;
            a = gradient.initial_a;
        }
        let before = (
            model.residual_calls.get(),
            model.load_calls.get(),
            model.load_vjp_calls.get(),
        );
        let got = recorded
            .pullback(
                (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
                &DIRECT,
                &Identity,
                budget(steps),
                &mut || false,
            )
            .unwrap();
        assert_eq!(got.initial_q, q);
        assert_eq!(got.initial_v, v);
        assert_eq!(got.initial_a, a);
        assert_eq!(got.parameters, parameters);
        let work = replay_costs(&costs);
        assert_eq!(got.replayed_steps, work.0);
        assert_eq!(
            model.residual_calls.get() - before.0,
            work.1,
            "duplicate primal leaf solve"
        );
        assert_eq!(model.load_calls.get() - before.1, work.0);
        assert_eq!(model.load_vjp_calls.get() - before.2, steps);
        assert_eq!(got.peak_checkpoints, budget(steps).checkpoints);
        assert_eq!(recorded.required_checkpoints(), budget(steps).checkpoints);
        for i in 0..2 {
            let delta = 1e-5 * (1.0 + model.p[i].abs());
            let mut plus = model.p;
            plus[i] += delta;
            let mut minus = model.p;
            minus[i] -= delta;
            let fd = (objective(plus, steps) - objective(minus, steps)) / (2.0 * delta);
            assert!(
                (got.parameters[i] - fd).abs() < 3e-7 * (1.0 + fd.abs()),
                "steps={steps}, parameter={i}, adjoint={}, FD={fd}",
                got.parameters[i]
            );
        }
    }
}

#[test]
fn forward_caps_and_cancellation_resume_the_same_accepted_prefix() {
    let model = Model::new([1.3, 0.6]);
    let mut initial = initial();
    let forcing = load(&model, initial.t);
    method().step(&mut initial, &model, &forcing).unwrap();
    let original = initial.clone();
    let mut recorded = RecordedStructural::new(method(), &model, &initial, config(7)).unwrap();
    assert_eq!(recorded.state().history, Vec::new());
    assert_eq!(initial, original);
    let report = recorded.advance(0, 7, &mut || false).unwrap();
    assert_eq!(report.status, StructuralRecordingStatus::StepLimit);
    assert_eq!(report.advanced, 0);
    let report = recorded.advance(7, 2, &mut || false).unwrap();
    assert_eq!(report.status, StructuralRecordingStatus::RecordLimit);
    assert_eq!(report.advanced, 2);
    assert_eq!(recorded.accepted_steps(), 2);
    assert!(matches!(
        recorded.pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(7),
            &mut || false
        ),
        Err(StructuralTrajectoryError::Incomplete)
    ));
    let before = recorded.state().clone();
    let residuals = model.residual_calls.get();
    let report = recorded
        .advance(7, 7, &mut || model.residual_calls.get() > residuals)
        .unwrap();
    assert_eq!(report.status, StructuralRecordingStatus::Cancelled);
    assert_eq!(report.advanced, 0);
    same_state(recorded.state(), &before);
    assert_eq!(recorded.accepted_steps(), 2);
    let report = recorded.advance(1, 7, &mut || false).unwrap();
    assert_eq!(report.status, StructuralRecordingStatus::StepLimit);
    assert_eq!(report.advanced, 1);
    let mut fork = recorded.clone();
    assert_eq!(recorded.advance(7, 7, &mut || false).unwrap().advanced, 4);
    while fork.advance(1, 7, &mut || false).unwrap().status != StructuralRecordingStatus::ReachedEnd
    {
    }
    same_state(recorded.state(), fork.state());
    assert_eq!(recorded.state().history, Vec::new());
    for _ in 0..7 {
        let forcing = load(&model, initial.t);
        method().step(&mut initial, &model, &forcing).unwrap();
    }
    same_state(recorded.state(), &initial);
}

#[test]
fn checkpoint_replay_and_cancellation_limits_are_atomic_and_retryable() {
    let model = Model::new([1.3, 0.6]);
    let mut recorded = recording(&model, 7);
    recorded.advance(7, 7, &mut || false).unwrap();
    let before = recorded.state().clone();
    let expected = recorded
        .pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(7),
            &mut || false,
        )
        .unwrap();
    let calls = model.residual_calls.get();
    let limited = StructuralReplayBudget {
        checkpoints: recorded.required_checkpoints() - 1,
        ..budget(7)
    };
    assert!(matches!(
        recorded.pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            limited,
            &mut || false
        ),
        Err(StructuralTrajectoryError::CheckpointLimit { .. })
    ));
    assert_eq!(model.residual_calls.get(), calls);
    for allowed in [0, expected.replayed_steps - 1] {
        let limited = StructuralReplayBudget {
            forward_steps: allowed,
            ..budget(7)
        };
        assert!(matches!(
            recorded.pullback(
                (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
                &DIRECT,
                &Identity,
                limited,
                &mut || false
            ),
            Err(StructuralTrajectoryError::ReplayLimit)
        ));
        same_state(recorded.state(), &before);
    }
    let residuals = model.residual_calls.get();
    assert!(matches!(
        recorded.pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(7),
            &mut || model.residual_calls.get() > residuals
        ),
        Err(StructuralTrajectoryError::Step(
            SecondOrderAdjointError::Step(TimeSolveError::Cancelled)
        ))
    ));
    same_state(recorded.state(), &before);
    let retried = recorded
        .pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(7),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retried.initial_q, expected.initial_q);
    assert_eq!(retried.initial_v, expected.initial_v);
    assert_eq!(retried.initial_a, expected.initial_a);
    assert_eq!(retried.parameters, expected.parameters);
}

#[test]
fn changed_loading_is_detected_before_a_leaf_derivative_is_used() {
    let model = Model::new([1.3, 0.6]);
    let mut recorded = recording(&model, 7);
    recorded.advance(7, 7, &mut || false).unwrap();
    let before = recorded.state().clone();
    model.fault.set(Fault::ChangedLoad);
    assert!(matches!(
        recorded.pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(7),
            &mut || false
        ),
        Err(StructuralTrajectoryError::ReplayMismatch { step: 0 })
    ));
    assert_eq!(model.load_vjp_calls.get(), 0);
    same_state(recorded.state(), &before);
}

#[test]
fn load_and_derivative_faults_refuse_without_mutating_the_recording() {
    let model = Model::new([1.3, 0.6]);
    let mut recorded = recording(&model, 3);
    recorded.advance(1, 3, &mut || false).unwrap();
    let before = recorded.state().clone();
    for fault in [Fault::LoadError, Fault::LoadIncomplete] {
        model.fault.set(fault);
        assert_fault(&recorded.advance(3, 3, &mut || false).unwrap_err(), fault);
        same_state(recorded.state(), &before);
        assert_eq!(recorded.accepted_steps(), 1);
    }
    model.fault.set(Fault::None);
    assert_eq!(
        recorded.advance(3, 3, &mut || false).unwrap().status,
        StructuralRecordingStatus::ReachedEnd
    );
    let before = recorded.state().clone();
    for fault in [
        Fault::LoadError,
        Fault::LoadIncomplete,
        Fault::DerivativeError,
        Fault::DerivativeIncomplete,
    ] {
        model.fault.set(fault);
        assert_fault(
            &recorded
                .pullback(
                    (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
                    &DIRECT,
                    &Identity,
                    budget(3),
                    &mut || false,
                )
                .unwrap_err(),
            fault,
        );
        same_state(recorded.state(), &before);
    }
    model.fault.set(Fault::None);
    recorded
        .pullback(
            (&TERMINAL.0, &TERMINAL.1, &TERMINAL.2),
            &DIRECT,
            &Identity,
            budget(3),
            &mut || false,
        )
        .unwrap();
}
