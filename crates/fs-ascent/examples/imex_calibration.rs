//! Calibrate a stiff two-body thermal network from synthetic sensor readings.
//!
//! Run `cargo run -p fs-ascent --features transient-design --example imex_calibration`.
//! Temperatures are excess above a fixed ambient, in kelvin. Heat capacities
//! are 0.01 and 2 J/K; ambient conductances are 0.1 and 0.2 W/K. The unknowns
//! are inter-body conductance (W/K), heater power (W), common sensor offset
//! (K), and the initial hot-body temperature (K). Unequal capacities make the
//! temperature operator nonsymmetric. Its fast mode needs implicit treatment
//! at the chosen 0.05 s step; the constant heater source is explicit.
//!
//! The existing dense IMEX solver generates noiseless readings. Calibration
//! uses matrix-free primal/transposed solves and checkpointed adjoints. This
//! demonstrates numerical parameter recovery, not experimental validation or
//! uncertainty/identifiability certification for a real thermal assembly.

use fs_ascent::SqpStop;
use fs_ascent::transient::imex::{
    ImexTransientConfig, ImexTransientFamily, ImexTransientModel, ImexTransientStudy,
};
use fs_solver::LinearOp;
use fs_time::stiff::adjoint::ImexVjp;
use fs_time::stiff::adjoint::trajectory::{
    ImexRecordingConfig, ImexReplayBudget, samples::SampleObjective,
};
use fs_time::stiff::{IdentityPreconditioner, Imex2, ImexSolveConfig, imex2_step};
use std::sync::Arc;

/// Synthetic conductance, heater power, sensor offset, and initial temperature.
pub const TRUTH: [f64; 4] = [1.2, 1.8, 0.15, 0.8];
/// Initial calibration point in the same physical coordinates as `TRUTH`.
pub const START: [f64; 4] = [0.7, 1.1, -0.2, 1.4];
/// Fixed sampling mesh spacing, in seconds.
pub const STEP_SECONDS: f64 = 0.05;
/// Number of steps spanning the 3.2-second observation interval.
pub const STEPS: usize = 64;
const CAPACITY: [f64; 2] = [0.01, 2.0];
const AMBIENT_CONDUCTANCE: [f64; 2] = [0.1, 0.2];
const BOUNDS: [[f64; 2]; 4] = [[0.2, 2.5], [0.1, 3.0], [-0.5, 0.5], [0.1, 2.0]];

/// Both temperature sensors on a shared, fixed endpoint timetable.
pub struct HeatReadings {
    /// Repeated endpoints select hot then cold sensors in declaration order.
    pub indices: Vec<usize>,
    targets: Arc<Vec<f64>>,
}

impl HeatReadings {
    /// Generate noiseless readings using the production dense IMEX solver.
    #[must_use]
    pub fn synthetic() -> Self {
        let indices: Vec<_> = [0, 1, 2, 4, 8, 16, 32, STEPS]
            .into_iter()
            .flat_map(|index| [index, index])
            .collect();
        let targets = Arc::new(dense_readings(&TRUTH, &indices));
        Self { indices, targets }
    }

    /// Independent forward value reference using the existing dense LU path.
    #[must_use]
    pub fn dense_loss(&self, point: &[f64; 4]) -> f64 {
        dense_readings(point, &self.indices)
            .iter()
            .zip(self.targets.iter())
            .map(|(value, target)| 0.5 * (value - target).powi(2))
            .sum()
    }
}

fn matrix(point: &[f64; 4]) -> [f64; 4] {
    let g = point[0];
    [
        -(g + AMBIENT_CONDUCTANCE[0]) / CAPACITY[0],
        g / CAPACITY[0],
        g / CAPACITY[1],
        -(g + AMBIENT_CONDUCTANCE[1]) / CAPACITY[1],
    ]
}

fn dense_readings(point: &[f64; 4], indices: &[usize]) -> Vec<f64> {
    let method = Imex2::new(&matrix(point), 2, STEP_SECONDS);
    let mut state = [point[3], 0.0];
    let mut step = 0;
    indices
        .iter()
        .enumerate()
        .map(|(sensor, &index)| {
            while step < index {
                imex2_step(&method, &mut state, &|_, out| {
                    out[0] = point[1] / CAPACITY[0];
                    out[1] = 0.0;
                });
                step += 1;
            }
            state[sensor % 2] + point[2]
        })
        .collect()
}

/// Immutable thermal operator, source, initial state, and sensor objective.
pub struct HeatModel {
    point: [f64; 4],
    initial: [f64; 2],
    targets: Arc<Vec<f64>>,
}

impl LinearOp for HeatModel {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, state: &[f64], out: &mut [f64]) {
        let a = matrix(&self.point);
        for i in 0..2 {
            out[i] = a[2 * i + 1].mul_add(state[1], a[2 * i] * state[0]);
        }
    }
    fn apply_transpose(&self, seed: &[f64], out: &mut [f64]) {
        let a = matrix(&self.point);
        for i in 0..2 {
            out[i] = a[2 + i].mul_add(seed[1], a[i] * seed[0]);
        }
    }
}

impl ImexVjp for HeatModel {
    fn parameter_count(&self) -> usize {
        4
    }
    fn nonlinear(&self, _: &[f64], out: &mut [f64]) {
        out[0] = self.point[1] / CAPACITY[0];
        out[1] = 0.0;
    }
    fn nonlinear_vjp(
        &self,
        _: &[f64],
        seed: &[f64],
        state_bar: &mut [f64],
        parameter_bar: &mut [f64],
    ) -> Result<(), String> {
        state_bar.fill(0.0);
        parameter_bar.fill(0.0);
        parameter_bar[1] = seed[0] / CAPACITY[0];
        Ok(())
    }
    fn linear_parameter_vjp(
        &self,
        state: &[f64],
        seed: &[f64],
        parameter_bar: &mut [f64],
    ) -> Result<(), String> {
        parameter_bar.fill(0.0);
        parameter_bar[0] = (state[0] - state[1]) * (-seed[0] / CAPACITY[0] + seed[1] / CAPACITY[1]);
        Ok(())
    }
}

impl SampleObjective for HeatModel {
    fn evaluate(
        &self,
        sample: usize,
        _: f64,
        state: &[f64],
        state_bar: &mut [f64],
        parameter_bar: &mut [f64],
    ) -> Result<f64, String> {
        let sensor = sample % 2;
        let target = self.targets.get(sample).ok_or("missing sensor reading")?;
        let residual = state[sensor] + self.point[2] - target;
        state_bar.fill(0.0);
        state_bar[sensor] = residual;
        parameter_bar.fill(0.0);
        parameter_bar[2] = residual;
        Ok(0.5 * residual * residual)
    }
}

impl ImexTransientModel for HeatModel {
    fn initial_values(&self) -> &[f64] {
        &self.initial
    }
    fn initial_vjp(&self, initial_bar: &[f64], parameter_bar: &mut [f64]) -> Result<(), String> {
        parameter_bar.fill(0.0);
        parameter_bar[3] = initial_bar[0];
        Ok(())
    }
}

impl ImexTransientFamily for HeatReadings {
    type Model = HeatModel;
    fn bounds(&self) -> &[[f64; 2]] {
        &BOUNDS
    }
    fn sample_indices(&self) -> &[usize] {
        &self.indices
    }
    fn instantiate(&self, point: &[f64]) -> Result<HeatModel, String> {
        let point: [f64; 4] = point
            .try_into()
            .map_err(|_| "expected four thermal parameters")?;
        Ok(HeatModel {
            point,
            initial: [point[3], 0.0],
            targets: self.targets.clone(),
        })
    }
}

/// Explicit small-system solve, recording, replay, and optimization budgets.
#[must_use]
pub fn config() -> ImexTransientConfig {
    ImexTransientConfig {
        start: 0.0,
        step: STEP_SECONDS,
        solve: ImexSolveConfig {
            tolerance: 1e-12,
            restart: 2,
            max_cycles: 4,
        },
        recording: ImexRecordingConfig {
            steps: STEPS,
            max_workspace_components: 512,
        },
        max_state_components: 2,
        max_samples: 16,
        max_forward_steps: STEPS,
        max_records: STEPS,
        replay: ImexReplayBudget {
            checkpoints: 8,
            forward_steps: 1024,
        },
        max_kkt_dimension: 12,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let readings = HeatReadings::synthetic();
    let primal = IdentityPreconditioner;
    let adjoint = IdentityPreconditioner;
    let mut study =
        ImexTransientStudy::new(&readings, &START, config(), &primal, &adjoint, &mut || {
            false
        })?;
    let initial = study.accepted().value;
    let report = study.run(1e-7, 100, 1000, &mut || false)?;
    let point = study.optimizer().point();
    println!(
        "source=synthetic-dense-imex model=two-body-heat stop={:?} iterations={} evaluations={}",
        report.stop,
        study.optimizer().iterations(),
        study.optimizer().evaluations()
    );
    println!(
        "conductance_W_per_K={:.9} heater_W={:.9} sensor_offset_K={:.9} initial_hot_K={:.9}",
        point[0], point[1], point[2], point[3]
    );
    println!(
        "objective_before={initial:.12e} objective_after={:.12e} observations={} replayed_steps={} peak_checkpoints={}",
        study.accepted().value,
        study.accepted().observations,
        study.accepted().replayed_steps,
        study.accepted().peak_checkpoints
    );
    println!("scope=numerical-parameter-recovery; no experimental-validation or uncertainty claim");
    if report.stop != SqpStop::Converged {
        return Err("calibration stopped before satisfying its local KKT tolerance".into());
    }
    Ok(())
}
