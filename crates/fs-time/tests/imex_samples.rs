//! G3/G4: sampled IMEX objectives against dense solves and full-storage reverse.
use fs_solver::LinearOp;
use fs_time::stiff::{
    IdentityPreconditioner, Imex2, ImexSolveConfig, ImexSolveError, ImexState, OperatorImex2,
    adjoint::{
        ImexAdjointError, ImexVjp,
        trajectory::{
            ImexRecordingConfig, ImexRecordingStatus, ImexReplayBudget, ImexTrajectoryError,
            RecordedImex2, samples::SampleObjective,
        },
    },
    imex2_step,
};
use std::cell::RefCell;

const START: f64 = 1.7;
const H: f64 = 0.08;

struct Model([f64; 2]);
impl Model {
    fn matrix(&self) -> [f64; 4] {
        [-self.0[0], 1.3, 0.4, -self.0[0] - 3.0]
    }
}
impl LinearOp for Model {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, u: &[f64], out: &mut [f64]) {
        let a = self.matrix();
        out[0] = a[0] * u[0] + a[1] * u[1];
        out[1] = a[2] * u[0] + a[3] * u[1];
    }
    fn apply_transpose(&self, u: &[f64], out: &mut [f64]) {
        let a = self.matrix();
        out[0] = a[0] * u[0] + a[2] * u[1];
        out[1] = a[1] * u[0] + a[3] * u[1];
    }
}
impl ImexVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }
    fn nonlinear(&self, u: &[f64], out: &mut [f64]) {
        out[0] = self.0[1] * u[0] * u[0] + 0.1 * u[1];
        out[1] = -0.2 * u[0] * u[1] + self.0[1];
    }
    fn nonlinear_vjp(
        &self,
        u: &[f64],
        seed: &[f64],
        ub: &mut [f64],
        pb: &mut [f64],
    ) -> Result<(), String> {
        ub[0] = 2.0 * self.0[1] * u[0] * seed[0] - 0.2 * u[1] * seed[1];
        ub[1] = 0.1 * seed[0] - 0.2 * u[0] * seed[1];
        pb.copy_from_slice(&[0.0, u[0] * u[0] * seed[0] + seed[1]]);
        Ok(())
    }
    fn linear_parameter_vjp(&self, u: &[f64], seed: &[f64], pb: &mut [f64]) -> Result<(), String> {
        pb.copy_from_slice(&[-u[0] * seed[0] - u[1] * seed[1], 0.0]);
        Ok(())
    }
}

fn method() -> OperatorImex2 {
    OperatorImex2::new(
        2,
        H,
        ImexSolveConfig {
            tolerance: 1e-13,
            restart: 2,
            max_cycles: 4,
        },
    )
}

fn initial(input: [f64; 4]) -> [f64; 2] {
    [input[0] + 0.1 * input[2], input[1] - 0.2 * input[3]]
}

fn recording<'a>(
    model: &'a Model,
    u: &[f64],
    steps: usize,
) -> RecordedImex2<'a, Model, IdentityPreconditioner> {
    RecordedImex2::new(
        method(),
        model,
        &IdentityPreconditioner,
        START,
        u,
        ImexRecordingConfig {
            steps,
            max_workspace_components: method().adjoint_workspace_components(2).unwrap(),
        },
    )
    .unwrap()
}

fn times(steps: usize) -> Vec<f64> {
    let mut result = vec![START];
    for step in 0..steps {
        result.push(result[step] + H);
    }
    result
}

fn term(sample: usize, time: f64, u: &[f64], p: [f64; 2]) -> (f64, [f64; 2], [f64; 2]) {
    let sensor = -0.4 + 0.13 * sample as f64;
    let weight = 0.5 + 0.1 * sample as f64;
    let target = 0.15 + 0.07 * sample as f64 + 0.03 * time + 0.1 * p[1];
    let error = u[0] + sensor * u[1] - target;
    let penalty = 0.025 * (sample + 1) as f64;
    let value = 0.5 * weight * error * error
        + penalty * p[0] * p[0]
        + 0.05 * p[0] * p[1]
        + 0.02 * time * p[1];
    (
        value,
        [weight * error, weight * error * sensor],
        [
            2.0 * penalty * p[0] + 0.05 * p[1],
            -0.1 * weight * error + 0.05 * p[0] + 0.02 * time,
        ],
    )
}

struct Objective {
    p: [f64; 2],
    expected_times: Vec<f64>,
    calls: RefCell<Vec<usize>>,
}
impl Objective {
    fn new(p: [f64; 2], indices: &[usize], steps: usize) -> Self {
        let endpoint_times = times(steps);
        Self {
            p,
            expected_times: indices.iter().map(|i| endpoint_times[*i]).collect(),
            calls: RefCell::new(Vec::new()),
        }
    }
}
impl SampleObjective for Objective {
    fn evaluate(
        &self,
        sample: usize,
        time: f64,
        u: &[f64],
        ub: &mut [f64],
        pb: &mut [f64],
    ) -> Result<f64, String> {
        assert_eq!(time.to_bits(), self.expected_times[sample].to_bits());
        self.calls.borrow_mut().push(sample);
        let (value, state, parameters) = term(sample, time, u, self.p);
        ub.copy_from_slice(&state);
        pb.copy_from_slice(&parameters);
        Ok(value)
    }
}

fn dense_value(input: [f64; 4], indices: &[usize], steps: usize) -> f64 {
    let model = Model([input[2], input[3]]);
    let dense = Imex2::new(&model.matrix(), 2, H);
    let mut states = vec![initial(input).to_vec()];
    for step in 0..steps {
        let mut next = states[step].clone();
        imex2_step(&dense, &mut next, &|u, out| model.nonlinear(u, out));
        states.push(next);
    }
    let endpoint_times = times(steps);
    indices
        .iter()
        .enumerate()
        .rev()
        .map(|(sample, endpoint)| {
            term(
                sample,
                endpoint_times[*endpoint],
                &states[*endpoint],
                model.0,
            )
            .0
        })
        .sum()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

#[test]
#[allow(clippy::too_many_lines)] // Keep the full-storage and dense oracles beside their comparisons.
fn samples_match_full_storage_and_dense_five_point_differences() {
    let input = [0.9, -0.4, 4.0, 0.7];
    let model = Model([input[2], input[3]]);
    let indices = [0, 0, 1, 3, 3, 5, 8, 8];
    let steps = 8;
    let mut tape = recording(&model, &initial(input), steps);
    assert_eq!(
        tape.advance(steps, steps, &mut || false).unwrap().status,
        ImexRecordingStatus::ReachedEnd
    );
    let objective = Objective::new(model.0, &indices, steps);
    let budget = ImexReplayBudget {
        checkpoints: tape.required_checkpoints(),
        forward_steps: 1000,
    };
    let actual = tape
        .pullback_samples(
            &indices,
            indices.len(),
            &objective,
            &IdentityPreconditioner,
            budget,
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.observations, indices.len());
    assert_eq!(
        *objective.calls.borrow(),
        (0..indices.len()).rev().collect::<Vec<_>>()
    );

    // Independent full-storage schedule: inject every endpoint's terms, then
    // reverse that step. Initial observations are added after the final VJP.
    let mut state = ImexState::new(START, &initial(input));
    let mut states = vec![state.u.clone()];
    for _ in 0..steps {
        method()
            .step(&mut state, &model, &IdentityPreconditioner, &|u, out| {
                model.nonlinear(u, out);
            })
            .unwrap();
        states.push(state.u.clone());
    }
    assert_eq!(bits(tape.state()), bits(&state.u));
    let endpoint_times = times(steps);
    let (mut value, mut state_bar, mut parameter_bar) = (0.0, vec![0.0; 2], vec![0.0; 2]);
    for endpoint in (0..=steps).rev() {
        for (sample, index) in indices.iter().enumerate().rev() {
            if *index == endpoint {
                let (v, ub, pb) =
                    term(sample, endpoint_times[endpoint], &states[endpoint], model.0);
                value += v;
                for i in 0..2 {
                    state_bar[i] += ub[i];
                    parameter_bar[i] += pb[i];
                }
            }
        }
        if endpoint > 0 {
            let gradient = method()
                .step_vjp(
                    &states[endpoint - 1],
                    &model,
                    &IdentityPreconditioner,
                    &IdentityPreconditioner,
                    &state_bar,
                    method().adjoint_workspace_components(2).unwrap(),
                    &mut || false,
                )
                .unwrap();
            state_bar = gradient.initial;
            for i in 0..2 {
                parameter_bar[i] += gradient.parameters[i];
            }
        }
    }
    assert_eq!(actual.value.to_bits(), value.to_bits());
    assert_eq!(bits(&actual.gradient.initial), bits(&state_bar));
    assert_eq!(bits(&actual.gradient.parameters), bits(&parameter_bar));
    close(actual.value, dense_value(input, &indices, steps));
    // Sampling adds no primal replays beyond the ordinary terminal sweep.
    let terminal = tape
        .pullback(
            &[0.0; 2],
            &[0.0; 2],
            &IdentityPreconditioner,
            budget,
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.gradient.replayed_steps, terminal.replayed_steps);
    assert_eq!(actual.gradient.peak_checkpoints, terminal.peak_checkpoints);

    // Initial-state parameterization is explicitly chained by the caller.
    let derivatives = [
        actual.gradient.initial[0],
        actual.gradient.initial[1],
        actual.gradient.parameters[0] + 0.1 * actual.gradient.initial[0],
        actual.gradient.parameters[1] - 0.2 * actual.gradient.initial[1],
    ];
    for i in 0..4 {
        let h = 2e-4 * (1.0 + input[i].abs());
        let values = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
            let mut perturbed = input;
            perturbed[i] += offset * h;
            dense_value(perturbed, &indices, steps)
        });
        let difference = (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h);
        close(derivatives[i], difference);
        assert!(difference.abs() > 1e-4, "uninformative derivative {i}");
    }
}

#[test]
fn initial_only_zero_step_objective_needs_no_checkpoints_or_replay() {
    let model = Model([4.0, 0.7]);
    let u = [0.9, -0.4];
    let tape = recording(&model, &u, 0);
    let indices = [0, 0, 0];
    let objective = Objective::new(model.0, &indices, 0);
    let actual = tape
        .pullback_samples(
            &indices,
            indices.len(),
            &objective,
            &IdentityPreconditioner,
            ImexReplayBudget {
                checkpoints: 0,
                forward_steps: 0,
            },
            &mut || false,
        )
        .unwrap();
    let (mut value, mut ub, mut pb) = (0.0, [0.0; 2], [0.0; 2]);
    for sample in (0..indices.len()).rev() {
        let (v, state, parameters) = term(sample, START, &u, model.0);
        value += v;
        for i in 0..2 {
            ub[i] += state[i];
            pb[i] += parameters[i];
        }
    }
    assert_eq!(actual.value.to_bits(), value.to_bits());
    assert_eq!(bits(&actual.gradient.initial), bits(&ub));
    assert_eq!(bits(&actual.gradient.parameters), bits(&pb));
    assert_eq!(actual.gradient.replayed_steps, 0);
    assert_eq!(actual.gradient.peak_checkpoints, 0);
    assert_eq!(actual.observations, indices.len());
    assert_eq!(*objective.calls.borrow(), [2, 1, 0]);
}

#[test]
fn malformed_or_incomplete_samples_refuse_before_objective_callbacks() {
    let model = Model([4.0, 0.7]);
    let mut tape = recording(&model, &[0.9, -0.4], 3);
    let objective = Objective::new(model.0, &[0, 1, 3], 3);
    let budget = ImexReplayBudget {
        checkpoints: 2,
        forward_steps: 100,
    };
    assert!(matches!(
        tape.pullback_samples(
            &[0, 1, 3],
            3,
            &objective,
            &IdentityPreconditioner,
            budget,
            &mut || false
        ),
        Err(ImexTrajectoryError::Incomplete)
    ));
    tape.advance(3, 3, &mut || false).unwrap();
    for (indices, cap) in [
        (&[][..], 3),
        (&[0, 3, 2][..], 3),
        (&[4][..], 3),
        (&[0, 1][..], 1),
    ] {
        assert!(matches!(
            tape.pullback_samples(
                indices,
                cap,
                &objective,
                &IdentityPreconditioner,
                budget,
                &mut || false
            ),
            Err(ImexTrajectoryError::InvalidInput(_))
        ));
    }
    assert!(objective.calls.borrow().is_empty());
}

struct BrokenObjective(u8);
impl SampleObjective for BrokenObjective {
    fn evaluate(
        &self,
        _sample: usize,
        _time: f64,
        _u: &[f64],
        ub: &mut [f64],
        pb: &mut [f64],
    ) -> Result<f64, String> {
        if self.0 == 0 {
            return Err("sensor unavailable".into());
        }
        if self.0 == 2 {
            ub[0] = 0.0;
        } else {
            ub.fill(0.0);
        }
        if self.0 == 3 {
            pb[0] = 0.0;
        } else {
            pb.fill(0.0);
        }
        Ok(if self.0 == 1 { f64::NAN } else { 0.0 })
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One accepted recording is retried after each named refusal.
fn observation_failures_budgets_and_every_cancellation_boundary_are_retryable() {
    let model = Model([4.0, 0.7]);
    let mut tape = recording(&model, &[0.9, -0.4], 3);
    tape.advance(3, 3, &mut || false).unwrap();
    let before = (
        bits(tape.state()),
        tape.time().to_bits(),
        tape.accepted_steps(),
    );
    let indices = [0, 1, 1, 3];
    let budget = ImexReplayBudget {
        checkpoints: 2,
        forward_steps: 100,
    };
    for mode in 0..=3 {
        assert!(matches!(
            tape.pullback_samples(
                &indices,
                indices.len(),
                &BrokenObjective(mode),
                &IdentityPreconditioner,
                budget,
                &mut || false
            ),
            Err(ImexTrajectoryError::Observation(_))
        ));
    }
    let objective = Objective::new(model.0, &indices, 3);
    let mut polls = 0;
    let reference = tape
        .pullback_samples(
            &indices,
            indices.len(),
            &objective,
            &IdentityPreconditioner,
            budget,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(matches!(
        tape.pullback_samples(
            &indices,
            indices.len(),
            &objective,
            &IdentityPreconditioner,
            ImexReplayBudget {
                checkpoints: 1,
                ..budget
            },
            &mut || false
        ),
        Err(ImexTrajectoryError::CheckpointLimit { .. })
    ));
    assert!(matches!(
        tape.pullback_samples(
            &indices,
            indices.len(),
            &objective,
            &IdentityPreconditioner,
            ImexReplayBudget {
                forward_steps: reference.gradient.replayed_steps - 1,
                ..budget
            },
            &mut || false
        ),
        Err(ImexTrajectoryError::ReplayLimit)
    ));
    for stop in 1..=polls {
        let mut calls = 0;
        let objective = Objective::new(model.0, &indices, 3);
        assert!(
            matches!(
                tape.pullback_samples(
                    &indices,
                    indices.len(),
                    &objective,
                    &IdentityPreconditioner,
                    budget,
                    &mut || {
                        calls += 1;
                        calls == stop
                    }
                ),
                Err(ImexTrajectoryError::Step(ImexAdjointError::Step(
                    ImexSolveError::Cancelled
                )))
            ),
            "cancellation poll {stop}"
        );
    }
    assert_eq!(
        (
            bits(tape.state()),
            tape.time().to_bits(),
            tape.accepted_steps()
        ),
        before
    );
    let retry = tape
        .pullback_samples(
            &indices,
            indices.len(),
            &Objective::new(model.0, &indices, 3),
            &IdentityPreconditioner,
            budget,
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.value.to_bits(), reference.value.to_bits());
    assert_eq!(retry.gradient, reference.gradient);
    assert_eq!(retry.observations, reference.observations);
}
