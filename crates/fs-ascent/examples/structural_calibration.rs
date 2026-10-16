//! Calibrate a forced two-body spring/damper model from q/v/a observations.
//!
//! Run `cargo run -p fs-ascent --features transient-design --example structural_calibration`.
//! Masses are 1 and 1.7 kg, with ground springs 0.8 and 1.1 N/m and ground
//! dampers 0.08 and 0.12 N s/m. Unknowns are coupling stiffness (N/m), coupling
//! damping (N s/m), forcing amplitude (N), initial displacement (m), and a
//! common displacement-sensor offset (m). The load varies with time on both
//! bodies; initial acceleration satisfies the physical equilibrium equation.
//!
//! Noiseless readings come from the existing dense generalized-alpha solver.
//! The fit uses matrix-free Newton steps and sampled checkpointed adjoints.
//! Residuals are nondimensionalized by 1 m, 1 m/s, and 1 m/s² for q/v/a.
//! This is numerical recovery for a lumped model, with no plate/shell,
//! experimentally validated damping, or apparatus-authority claim.

use fs_ascent::SqpStop;
use fs_ascent::transient::structural::{
    StructuralTransientConfig, StructuralTransientFamily, StructuralTransientModel,
    StructuralTransientStudy,
};
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::initialization::{
    InitialSolveConfig, second_order_acceleration, second_order_acceleration_vjp,
};
use fs_time::galpha::second_order_adjoint::trajectory::{
    StructuralRecordingConfig, StructuralReplayBudget, StructuralTrajectoryModel,
    samples::StructuralSampleObjective,
};
use fs_time::galpha::second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderVjp};
use fs_time::galpha::{
    GeneralizedAlpha, ImplicitSolveConfig, SecondOrderProblem, SecondOrderState, galpha_step,
};
use fs_time::stiff::IdentityPreconditioner;
use std::sync::Arc;

/// Synthetic coupling stiffness, damping, load, initial displacement and offset.
pub const TRUTH: [f64; 5] = [1.2, 0.18, 1.1, 0.45, 0.04];
/// Starting values in the same physical coordinates as `TRUTH`.
pub const START: [f64; 5] = [0.85, 0.30, 0.8, 0.25, -0.06];
/// Fixed time step, in seconds.
pub const STEP_SECONDS: f64 = 0.08;
/// Number of steps across the 3.84-second observation interval.
pub const STEPS: usize = 48;
const RHO_INF: f64 = 0.65;
const MASS: [f64; 2] = [1.0, 1.7];
const GROUND_K: [f64; 2] = [0.8, 1.1];
const GROUND_C: [f64; 2] = [0.08, 0.12];
const INITIAL_V: [f64; 2] = [0.2, -0.1];
const INITIAL_SOLVE: InitialSolveConfig = InitialSolveConfig {
    restart: 2,
    max_cycles: 4,
    tolerance: 1e-12,
};
const INITIAL_WORKSPACE: usize = 256;
const BOUNDS: [[f64; 2]; 5] = [
    [0.5, 2.5],
    [0.04, 0.65],
    [0.3, 2.0],
    [0.1, 0.8],
    [-0.25, 0.25],
];

fn coupled(ground: [f64; 2], coupling: f64) -> [f64; 4] {
    [
        ground[0] + coupling,
        -coupling,
        -coupling,
        ground[1] + coupling,
    ]
}

fn multiply(matrix: [f64; 4], x: &[f64], out: &mut [f64]) {
    for i in 0..2 {
        out[i] = matrix[2 * i + 1].mul_add(x[1], matrix[2 * i] * x[0]);
    }
}

fn load_shape(t: f64) -> [f64; 2] {
    [1.0 + 0.4 * (1.7 * t).sin(), 0.25 * (0.9 * t).cos()]
}

fn dense_initial(point: &[f64; 5]) -> SecondOrderState {
    let q = [point[3], -0.15];
    let (mut damping, mut stiffness) = ([0.0; 2], [0.0; 2]);
    multiply(coupled(GROUND_C, point[1]), &INITIAL_V, &mut damping);
    multiply(coupled(GROUND_K, point[0]), &q, &mut stiffness);
    let force = load_shape(0.0);
    let a = std::array::from_fn::<_, 2, _>(|i| {
        (point[2] * force[i] - damping[i] - stiffness[i]) / MASS[i]
    });
    SecondOrderState::new(0.0, &q, &INITIAL_V, &a)
}

fn reading(sample: usize, state: &SecondOrderState, point: &[f64; 5]) -> f64 {
    let i = sample % 2;
    match sample % 6 {
        0 | 1 => state.q[i] + point[4],
        2 | 3 => state.v[i],
        _ => state.a[i],
    }
}

/// Displacement, velocity and acceleration on one fixed endpoint timetable.
pub struct StructuralReadings {
    /// Six repeated endpoints select q0/q1/v0/v1/a0/a1 in declaration order.
    pub indices: Vec<usize>,
    targets: Arc<Vec<f64>>,
}

impl StructuralReadings {
    /// Generate noiseless readings with the independent dense stepping path.
    #[must_use]
    pub fn synthetic() -> Self {
        let indices: Vec<_> = [0, 2, 4, 8, 16, 24, 36, STEPS]
            .into_iter()
            .flat_map(|i| [i; 6])
            .collect();
        let targets = Arc::new(dense_readings(&TRUTH, &indices));
        Self { indices, targets }
    }

    /// Independent forward reference for the complete least-squares objective.
    #[must_use]
    pub fn dense_loss(&self, point: &[f64; 5]) -> f64 {
        dense_readings(point, &self.indices)
            .iter()
            .zip(self.targets.iter())
            .map(|(actual, target)| 0.5 * (actual - target).powi(2))
            .sum()
    }
}

fn dense_readings(point: &[f64; 5], indices: &[usize]) -> Vec<f64> {
    let method = GeneralizedAlpha::new(
        &[MASS[0], 0.0, 0.0, MASS[1]],
        &coupled(GROUND_C, point[1]),
        &coupled(GROUND_K, point[0]),
        2,
        STEP_SECONDS,
        RHO_INF,
    );
    let mut state = dense_initial(point);
    let alpha_f = RHO_INF / (1.0 + RHO_INF);
    indices
        .iter()
        .enumerate()
        .map(|(sample, &index)| {
            while state.steps < index {
                let stage_time = (1.0 - alpha_f).mul_add(STEP_SECONDS, state.t);
                let force = load_shape(stage_time).map(|value| point[2] * value);
                galpha_step(&method, &mut state.q, &mut state.v, &mut state.a, &force);
                state.t += STEP_SECONDS;
                state.steps += 1;
            }
            reading(sample, &state, point)
        })
        .collect()
}

/// Immutable two-body residual, loads, consistent initial state, and sensors.
pub struct StructuralModel {
    point: [f64; 5],
    initial: SecondOrderState,
    targets: Arc<Vec<f64>>,
}

impl SecondOrderProblem for StructuralModel {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        for i in 0..2 {
            output[i] = MASS[i] * input[i];
        }
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        multiply(coupled(GROUND_C, self.point[1]), input, output);
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        multiply(coupled(GROUND_K, self.point[0]), q, output);
    }
    fn tangent_apply(&self, _: &[f64], direction: &[f64], output: &mut [f64]) {
        self.internal_force(direction, output);
    }
}

impl SecondOrderVjp for StructuralModel {
    fn parameter_count(&self) -> usize {
        5
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.mass_apply(seed, output);
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.damping_apply(seed, output);
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.tangent_apply(q, seed, output);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        _: &[f64],
        seed: &[f64],
        pbar: &mut [f64],
    ) -> Result<(), String> {
        pbar.fill(0.0);
        pbar[0] = (q[0] - q[1]) * (seed[0] - seed[1]);
        pbar[1] = (v[0] - v[1]) * (seed[0] - seed[1]);
        Ok(())
    }
}

impl StructuralTrajectoryModel for StructuralModel {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        for (out, value) in output.iter_mut().zip(load_shape(time)) {
            *out = self.point[2] * value;
        }
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], pbar: &mut [f64]) -> Result<(), String> {
        let shape = load_shape(time);
        pbar.fill(0.0);
        pbar[2] = shape[0] * seed[0] + shape[1] * seed[1];
        Ok(())
    }
}

impl StructuralSampleObjective for StructuralModel {
    fn evaluate(
        &self,
        sample: usize,
        state: &SecondOrderState,
        state_bar: (&mut [f64], &mut [f64], &mut [f64]),
        pbar: &mut [f64],
    ) -> Result<f64, String> {
        let target = self
            .targets
            .get(sample)
            .ok_or("missing structural reading")?;
        let residual = reading(sample, state, &self.point) - target;
        let (qbar, vbar, abar) = state_bar;
        qbar.fill(0.0);
        vbar.fill(0.0);
        abar.fill(0.0);
        pbar.fill(0.0);
        let i = sample % 2;
        match sample % 6 {
            0 | 1 => {
                qbar[i] = residual;
                pbar[4] = residual;
            }
            2 | 3 => vbar[i] = residual,
            _ => abar[i] = residual,
        }
        Ok(0.5 * residual * residual)
    }
}

impl StructuralTransientModel for StructuralModel {
    fn initial_state(&self) -> &SecondOrderState {
        &self.initial
    }
    fn initial_vjp(
        &self,
        qbar: &[f64],
        _: &[f64],
        abar: &[f64],
        pbar: &mut [f64],
    ) -> Result<(), String> {
        let force = load_shape(0.0).map(|value| self.point[2] * value);
        // The consumer polls around this bounded callback. Its trait does not
        // carry a cancellation token into the small initialization solve.
        let pullback = second_order_acceleration_vjp(
            self,
            &self.initial.q,
            &self.initial.v,
            &force,
            abar,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            INITIAL_SOLVE,
            INITIAL_WORKSPACE,
            &mut || false,
        )
        .map_err(|error| error.to_string())?;
        self.forcing_vjp(0.0, &pullback.forcing, pbar)?;
        for (total, partial) in pbar.iter_mut().zip(pullback.parameters) {
            *total += partial;
        }
        pbar[3] += qbar[0] + pullback.initial_q[0];
        // Initial velocity and the second displacement are parameter-independent.
        Ok(())
    }
}

impl StructuralTransientFamily for StructuralReadings {
    type Model = StructuralModel;
    fn bounds(&self) -> &[[f64; 2]] {
        &BOUNDS
    }
    fn sample_indices(&self) -> &[usize] {
        &self.indices
    }
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String> {
        let point = point
            .try_into()
            .map_err(|_| "expected five structural parameters")?;
        let mut model = StructuralModel {
            point,
            initial: SecondOrderState::new(0.0, &[point[3], -0.15], &INITIAL_V, &[0.0; 2]),
            targets: self.targets.clone(),
        };
        let force = load_shape(0.0).map(|value| point[2] * value);
        model.initial.a = second_order_acceleration(
            &model,
            &model.initial.q,
            &model.initial.v,
            &force,
            &IdentityPreconditioner,
            INITIAL_SOLVE,
            INITIAL_WORKSPACE,
            &mut || false,
        )
        .map_err(|error| error.to_string())?
        .value;
        Ok(model)
    }
}

/// Explicit solve, observation, recording, replay and dense optimization limits.
#[must_use]
pub fn config() -> StructuralTransientConfig {
    StructuralTransientConfig {
        start: 0.0,
        step: STEP_SECONDS,
        rho_inf: RHO_INF,
        solve: ImplicitSolveConfig {
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
        recording: StructuralRecordingConfig {
            steps: STEPS,
            adjoint: SecondOrderAdjointConfig {
                restart: 2,
                max_cycles: 4,
                tolerance: 1e-12,
            },
            max_workspace_components: 1024,
        },
        max_state_components: 2,
        max_samples: 48,
        max_forward_steps: STEPS,
        max_records: STEPS,
        replay: StructuralReplayBudget {
            checkpoints: 8,
            forward_steps: 1024,
        },
        max_kkt_dimension: 15,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let data = StructuralReadings::synthetic();
    let adjoint = IdentityPreconditioner;
    let mut study =
        StructuralTransientStudy::new(&data, &START, config(), &adjoint, &mut || false)?;
    let initial = study.accepted().value;
    let report = study.run(1e-7, 100, 1000, &mut || false)?;
    let p = study.optimizer().point();
    println!(
        "source=synthetic-dense-generalized-alpha model=two-body-structural stop={:?} iterations={} evaluations={}",
        report.stop,
        study.optimizer().iterations(),
        study.optimizer().evaluations()
    );
    println!(
        "coupling_N_per_m={:.9} damping_Ns_per_m={:.9} load_N={:.9} initial_q_m={:.9} sensor_offset_m={:.9}",
        p[0], p[1], p[2], p[3], p[4]
    );
    println!(
        "objective_before={initial:.12e} objective_after={:.12e} observations={} replayed_steps={} peak_checkpoints={}",
        study.accepted().value,
        study.accepted().observations,
        study.accepted().replayed_steps,
        study.accepted().peak_checkpoints
    );
    println!(
        "scope=numerical-lumped-model-recovery; no plate/shell, physical-damping-validation or apparatus-authority claim"
    );
    if report.stop != SqpStop::Converged {
        return Err("calibration stopped before its local KKT tolerance".into());
    }
    Ok(())
}
