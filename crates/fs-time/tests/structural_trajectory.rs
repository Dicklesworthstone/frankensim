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
            samples::StructuralSampleObjective,
        },
    },
};
use std::cell::{Cell, RefCell};

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

#[derive(Clone, Copy, PartialEq, Eq)]
enum SampleFault {
    None,
    Callback,
    AccelerationIncomplete,
    ParameterIncomplete,
}

struct Sensors {
    p: [f64; 2],
    calls: RefCell<Vec<(usize, usize)>>,
    fault: Cell<SampleFault>,
}

impl Sensors {
    fn new(p: [f64; 2]) -> Self {
        Self {
            p,
            calls: RefCell::new(Vec::new()),
            fault: Cell::new(SampleFault::None),
        }
    }
}

fn sample_value(sample: usize, state: &SecondOrderState, p: [f64; 2]) -> f64 {
    let weight = 1.0 + 0.1 * sample as f64;
    0.5 * weight * state.q[0].powi(2)
        + 0.2 * state.v[1].powi(2)
        + 0.03 * state.a[0].powi(2)
        + p[0] * state.q[1]
        + 0.1 * p[1] * state.a[1]
        + 0.01 * p[0].powi(2)
}

impl StructuralSampleObjective for Sensors {
    fn evaluate(
        &self,
        sample: usize,
        state: &SecondOrderState,
        state_bar: (&mut [f64], &mut [f64], &mut [f64]),
        parameter_bar: &mut [f64],
    ) -> Result<f64, String> {
        self.calls.borrow_mut().push((sample, state.steps));
        // Sample four is reached after a valid later sample has accumulated.
        let fault = if sample == 4 {
            self.fault.get()
        } else {
            SampleFault::None
        };
        if fault == SampleFault::Callback {
            return Err("sensor unavailable".into());
        }
        state_bar.0[0] = (1.0 + 0.1 * sample as f64) * state.q[0];
        state_bar.0[1] = self.p[0];
        state_bar.1[0] = 0.0;
        state_bar.1[1] = 0.4 * state.v[1];
        state_bar.2[0] = 0.06 * state.a[0];
        if fault != SampleFault::AccelerationIncomplete {
            state_bar.2[1] = 0.1 * self.p[1];
        }
        parameter_bar[0] = state.q[1] + 0.02 * self.p[0];
        if fault != SampleFault::ParameterIncomplete {
            parameter_bar[1] = 0.1 * state.a[1];
        }
        Ok(sample_value(sample, state, self.p))
    }
}

fn sampled_value(p: [f64; 2], indices: &[usize]) -> f64 {
    let (states, _) = full_storage(&Model::new(p), *indices.last().unwrap());
    indices
        .iter()
        .enumerate()
        .rev()
        .map(|(sample, endpoint)| sample_value(sample, &states[*endpoint], p))
        .sum()
}

#[test]
#[allow(clippy::too_many_lines)] // One full-storage reference and its independent FD check.
fn repeated_structural_sensors_match_full_storage_and_parameter_differences_in_one_sweep() {
    let model = Model::new([1.3, 0.6]);
    let indices = [0, 0, 1, 4, 4, 7];
    let mut recorded = recording(&model, 7);
    recorded.advance(7, 7, &mut || false).unwrap();
    let (states, costs) = full_storage(&model, 7);
    let (mut q, mut v, mut a) = (vec![0.0; 2], vec![0.0; 2], vec![0.0; 2]);
    let mut parameters = vec![0.0; 2];
    let mut value = 0.0;
    for endpoint in (0..=7).rev() {
        let state = &states[endpoint];
        for (sample, _) in indices
            .iter()
            .enumerate()
            .rev()
            .filter(|(_, index)| **index == endpoint)
        {
            value += sample_value(sample, state, model.p);
            q[0] += (1.0 + 0.1 * sample as f64) * state.q[0];
            q[1] += model.p[0];
            v[1] += 0.4 * state.v[1];
            a[0] += 0.06 * state.a[0];
            a[1] += 0.1 * model.p[1];
            parameters[0] += state.q[1] + 0.02 * model.p[0];
            parameters[1] += 0.1 * state.a[1];
        }
        if endpoint == 0 {
            break;
        }
        let previous = &states[endpoint - 1];
        let gradient = method()
            .step_vjp(
                previous,
                &model,
                &load(&model, previous.t),
                (&q, &v, &a),
                &Identity,
                config(7).adjoint,
                config(7).max_workspace_components,
                &mut || false,
            )
            .unwrap();
        let mut forcing = [f64::NAN; 2];
        model
            .forcing_vjp(
                method().forcing_time(previous.t).unwrap(),
                &gradient.forcing,
                &mut forcing,
            )
            .unwrap();
        for i in 0..2 {
            parameters[i] += gradient.parameters[i];
            parameters[i] += forcing[i];
        }
        q = gradient.initial_q;
        v = gradient.initial_v;
        a = gradient.initial_a;
    }
    let sensors = Sensors::new(model.p);
    let before = (model.residual_calls.get(), model.load_calls.get());
    let result = recorded
        .pullback_samples(
            &indices,
            indices.len(),
            &sensors,
            &Identity,
            budget(7),
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.value.to_bits(), value.to_bits());
    assert_eq!(result.gradient.initial_q, q);
    assert_eq!(result.gradient.initial_v, v);
    assert_eq!(result.gradient.initial_a, a);
    assert_eq!(result.gradient.parameters, parameters);
    assert_eq!(result.observations, indices.len());
    assert_eq!(
        *sensors.calls.borrow(),
        indices
            .iter()
            .enumerate()
            .rev()
            .map(|(sample, endpoint)| (sample, states[*endpoint].steps))
            .collect::<Vec<_>>()
    );
    let work = replay_costs(&costs);
    assert_eq!(result.gradient.replayed_steps, work.0);
    assert_eq!(result.gradient.peak_checkpoints, budget(7).checkpoints);
    assert_eq!(model.residual_calls.get() - before.0, work.1);
    assert_eq!(model.load_calls.get() - before.1, work.0);
    for i in 0..2 {
        let delta = 1e-5 * (1.0 + model.p[i].abs());
        let mut plus = model.p;
        plus[i] += delta;
        let mut minus = model.p;
        minus[i] -= delta;
        let fd = (sampled_value(plus, &indices) - sampled_value(minus, &indices)) / (2.0 * delta);
        assert!(
            (result.gradient.parameters[i] - fd).abs() < 3e-7 * (1.0 + fd.abs()),
            "parameter {i}: adjoint={}, FD={fd}",
            result.gradient.parameters[i]
        );
    }
}

#[test]
fn repeated_initial_sensors_need_no_primal_steps_or_checkpoints() {
    let model = Model::new([1.3, 0.6]);
    let recorded = recording(&model, 0);
    let sensors = Sensors::new(model.p);
    let result = recorded
        .pullback_samples(&[0, 0], 2, &sensors, &Identity, budget(0), &mut || false)
        .unwrap();
    let state = initial();
    assert_eq!(
        result.value,
        sample_value(1, &state, model.p) + sample_value(0, &state, model.p)
    );
    assert_eq!(result.observations, 2);
    assert_eq!(
        result.gradient.initial_q,
        [1.1 * state.q[0] + state.q[0], 2.0 * model.p[0]]
    );
    assert_eq!(result.gradient.initial_v, [0.0, 0.8 * state.v[1]]);
    assert_eq!(
        result.gradient.initial_a,
        [0.12 * state.a[0], 0.2 * model.p[1]]
    );
    assert_eq!(
        result.gradient.parameters,
        [2.0 * (state.q[1] + 0.02 * model.p[0]), 0.2 * state.a[1]]
    );
    assert_eq!(
        *sensors.calls.borrow(),
        [(1, state.steps), (0, state.steps)]
    );
    assert_eq!(result.gradient.replayed_steps, 0);
    assert_eq!(result.gradient.peak_checkpoints, 0);
    assert_eq!(model.residual_calls.get(), 0);
    assert_eq!(model.load_calls.get(), 0);
}

#[test]
fn sampled_objective_bounds_faults_and_cancellation_preserve_retryable_recording() {
    let model = Model::new([1.3, 0.6]);
    let indices = [0, 0, 1, 4, 4, 7];
    let mut recorded = recording(&model, 7);
    recorded.advance(7, 7, &mut || false).unwrap();
    let before = recorded.state().clone();
    let sensors = Sensors::new(model.p);
    let calls = model.residual_calls.get();
    for (invalid, cap) in [
        (&[][..], 0),
        (&[4, 1][..], 2),
        (&[0, 8][..], 2),
        (&indices[..], indices.len() - 1),
    ] {
        assert!(matches!(
            recorded.pullback_samples(invalid, cap, &sensors, &Identity, budget(7), &mut || false),
            Err(StructuralTrajectoryError::InvalidInput(_))
        ));
    }
    assert_eq!(*sensors.calls.borrow(), Vec::new());
    assert_eq!(model.residual_calls.get(), calls);
    let expected = recorded
        .pullback_samples(
            &indices,
            indices.len(),
            &sensors,
            &Identity,
            budget(7),
            &mut || false,
        )
        .unwrap();
    for fault in [
        SampleFault::Callback,
        SampleFault::AccelerationIncomplete,
        SampleFault::ParameterIncomplete,
    ] {
        sensors.fault.set(fault);
        sensors.calls.borrow_mut().clear();
        assert!(matches!(
            recorded.pullback_samples(
                &indices,
                indices.len(),
                &sensors,
                &Identity,
                budget(7),
                &mut || false
            ),
            Err(StructuralTrajectoryError::Observation(_))
        ));
        assert_eq!(
            *sensors.calls.borrow(),
            [(5, before.steps), (4, initial().steps + 4)]
        );
        same_state(recorded.state(), &before);
    }
    sensors.fault.set(SampleFault::None);
    sensors.calls.borrow_mut().clear();
    assert!(matches!(
        recorded.pullback_samples(
            &indices,
            indices.len(),
            &sensors,
            &Identity,
            budget(7),
            &mut || sensors.calls.borrow().len() >= 2
        ),
        Err(StructuralTrajectoryError::Step(
            SecondOrderAdjointError::Step(TimeSolveError::Cancelled)
        ))
    ));
    same_state(recorded.state(), &before);
    sensors.calls.borrow_mut().clear();
    let retried = recorded
        .pullback_samples(
            &indices,
            indices.len(),
            &sensors,
            &Identity,
            budget(7),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retried.value.to_bits(), expected.value.to_bits());
    assert_eq!(retried.gradient, expected.gradient);
    assert_eq!(retried.observations, indices.len());
}
