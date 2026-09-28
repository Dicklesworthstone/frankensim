//! Fit stiffness, moving-force amplitude, and a displacement-sensor offset.
//!
//! Run `cargo run -p fs-ascent --features transient-design --example plate_calibration`.
//! The shared 32-triangle DKT fixture retains all 27 free spatial DOFs. Its
//! mass and Rayleigh coefficients are fixed: freeing mass, stiffness and load
//! scales together would introduce a scale ambiguity. The offset coordinate is
//! measured in units of 100 micrometres; force amplitude is in newtons.
//!
//! Synthetic readings use dense generalized-alpha stepping of the assembled
//! pencil. Each SQP trial instead uses the sparse production plate adapter,
//! consistent initial acceleration, and one sampled checkpointed reverse sweep.
//! This demonstrates numerical recovery on a fixed mesh, not experimental
//! material identification or convergence to a continuum plate solution.

#[allow(dead_code)]
#[path = "../../fs-plate-transient/examples/moving_plate.rs"]
mod moving;

use fs_ascent::sqp::{SqpSample, SqpState, SqpStop};
use fs_plate_transient::{PlateDynamics, PlateDynamicsParameters, PlateLoad};
use fs_plate::{PlateMesh, PlateModel};
use fs_time::galpha::initialization::{second_order_acceleration, second_order_acceleration_vjp};
use fs_time::galpha::second_order_adjoint::trajectory::{
    RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
    StructuralReplayBudget, StructuralTrajectoryModel, samples::StructuralSampleObjective,
};
use fs_time::galpha::{GeneralizedAlpha, SecondOrderState, galpha_step};

/// Fixed synthetic stiffness scale, force in N, and offset in 100-micrometre units.
pub const TRUTH: [f64; 3] = [1.15, 1.1, 0.08];
/// Initial SQP point in the same coordinates as `TRUTH`.
pub const START: [f64; 3] = [0.95, 0.75, -0.1];
/// Displacement normalization and sensor-offset coordinate, in metres.
pub const DISPLACEMENT_UNIT: f64 = 1e-4;
const BOUNDS: [[f64; 2]; 3] = [[0.6, 1.7], [0.3, 2.0], [-0.5, 0.5]];
const INITIAL_WORKSPACE: usize = 16_384;

fn parameters(point: &[f64]) -> PlateDynamicsParameters {
    PlateDynamicsParameters {
        stiffness_scale: point[0],
        mass_scale: 1.0,
        mass_damping_per_s: 4.0,
        stiffness_damping_s: 1e-5,
    }
}

// The second declared load coordinate is sensor-only: its forcing derivative
// is exactly zero, while observations write its direct partial in slot five.
struct CalibrationLoad<'a>(moving::MovingLoad<'a>);
impl PlateLoad for CalibrationLoad<'_> {
    fn parameter_count(&self) -> usize {
        2
    }
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        self.0.forcing(time, output)
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        output.fill(0.0);
        self.0.forcing_vjp(time, seed, &mut output[..1])
    }
}

/// Fixed assembled spatial problem, observation schedule and synthetic data.
pub struct PlateCalibration {
    /// Spatial DKT mesh; held fixed during every parameter evaluation.
    pub mesh: PlateMesh,
    /// Actual assembled stiffness, lumped mass, and clamped support map.
    pub plate: PlateModel,
    /// P1 displacement probe weights in reduced spatial coordinates.
    pub sensor: Vec<f64>,
    /// Accepted endpoint indices; zero observes the initial state.
    pub indices: Vec<usize>,
    /// Synthetic sensor displacement in metres, in timetable order.
    pub targets: Vec<f64>,
}

struct Observations<'a> {
    data: &'a PlateCalibration,
    offset: f64,
}
impl StructuralSampleObjective for Observations<'_> {
    fn evaluate(
        &self,
        sample: usize,
        state: &SecondOrderState,
        state_bar: (&mut [f64], &mut [f64], &mut [f64]),
        parameter_bar: &mut [f64],
    ) -> Result<f64, String> {
        let target = self
            .data
            .targets
            .get(sample)
            .ok_or("missing plate observation")?;
        let measured = self
            .data
            .sensor
            .iter()
            .zip(&state.q)
            .map(|(w, q)| w * q)
            .sum::<f64>();
        let residual = (measured + DISPLACEMENT_UNIT * self.offset - target) / DISPLACEMENT_UNIT;
        let (qbar, vbar, abar) = state_bar;
        for (bar, weight) in qbar.iter_mut().zip(&self.data.sensor) {
            *bar = residual * weight / DISPLACEMENT_UNIT;
        }
        vbar.fill(0.0);
        abar.fill(0.0);
        parameter_bar.fill(0.0);
        parameter_bar[5] = residual;
        Ok(0.5 * residual * residual)
    }
}

impl PlateCalibration {
    /// Assemble a real clamped plate and generate dense-reference readings.
    pub fn synthetic() -> Result<Self, Box<dyn std::error::Error>> {
        let (mesh, plate) = moving::build_plate()?;
        let sensor = moving::sensor_weights(&mesh, &plate)?;
        let mut data = Self {
            mesh,
            plate,
            sensor,
            indices: vec![0, 1, 2, 4, 6, 8, 10, moving::STEPS],
            targets: Vec::new(),
        };
        data.targets = data.dense_values(&TRUTH)?;
        Ok(data)
    }

    fn load(&self, point: &[f64]) -> CalibrationLoad<'_> {
        CalibrationLoad(moving::MovingLoad {
            mesh: &self.mesh,
            plate: &self.plate,
            amplitude: point[1],
        })
    }

    /// Independent dense objective for finite differences of all three unknowns.
    pub fn dense_loss(&self, point: &[f64; 3]) -> Result<f64, String> {
        Ok(self
            .dense_values(point)?
            .iter()
            .zip(&self.targets)
            .map(|(value, target)| 0.5 * ((value - target) / DISPLACEMENT_UNIT).powi(2))
            .sum())
    }

    fn dense_values(&self, point: &[f64; 3]) -> Result<Vec<f64>, String> {
        let n = self.plate.free;
        let p = parameters(point);
        let (mut mass, mut damping, mut stiffness) =
            (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]);
        for i in 0..n {
            for j in 0..n {
                mass[i * n + j] = self.plate.m.get(i, j);
                stiffness[i * n + j] = p.stiffness_scale * self.plate.k.get(i, j);
                damping[i * n + j] = p.mass_damping_per_s * mass[i * n + j]
                    + p.stiffness_damping_s * stiffness[i * n + j];
            }
        }
        let method =
            GeneralizedAlpha::new(&mass, &damping, &stiffness, n, moving::STEP, moving::RHO);
        let load = self.load(point);
        let mut force = vec![0.0; n];
        load.forcing(0.0, &mut force)?;
        // Independent diagonal equilibrium solve, with parameter-independent q0=v0=0.
        let acceleration: Vec<_> = force
            .iter()
            .enumerate()
            .map(|(i, f)| f / mass[i * n + i])
            .collect();
        let mut state = SecondOrderState::new(0.0, &vec![0.0; n], &vec![0.0; n], &acceleration);
        let alpha_f = moving::RHO / (1.0 + moving::RHO);
        let mut values = Vec::with_capacity(self.indices.len());
        for &index in &self.indices {
            while state.steps < index {
                load.forcing((1.0 - alpha_f).mul_add(moving::STEP, state.t), &mut force)?;
                galpha_step(&method, &mut state.q, &mut state.v, &mut state.a, &force);
                state.t += moving::STEP;
                state.steps += 1;
            }
            values.push(
                self.sensor
                    .iter()
                    .zip(&state.q)
                    .map(|(w, q)| w * q)
                    .sum::<f64>()
                    + point[2] * DISPLACEMENT_UNIT,
            );
        }
        Ok(values)
    }

    /// One bounded sparse forward/reverse trial, suitable for existing `SqpState`.
    /// `None` denotes an out-of-box trial; numerical refusal remains an error.
    #[allow(clippy::too_many_lines)]
    pub fn evaluate<C: FnMut() -> bool>(
        &self,
        point: &[f64],
        cancelled: &mut C,
    ) -> Result<Option<SqpSample>, String> {
        if cancelled() {
            return Err("plate calibration cancelled".into());
        }
        if point.len() != 3 || point.iter().any(|x| !x.is_finite()) {
            return Err("three finite plate parameters required".into());
        }
        if point
            .iter()
            .zip(BOUNDS)
            .any(|(x, b)| *x < b[0] || *x > b[1])
        {
            return Ok(None);
        }
        let load = self.load(point);
        let mut budget = moving::dynamics_budget(&self.plate);
        budget.max_load_parameters = 2;
        let dynamics = PlateDynamics::new(&self.plate, parameters(point), &load, budget, cancelled)
            .map_err(|e| e.to_string())?;
        let n = self.plate.free;
        let zero = vec![0.0; n];
        let mut force = zero.clone();
        dynamics.forcing(0.0, &mut force)?;
        let mass_preconditioner = dynamics.mass_preconditioner();
        let initial = second_order_acceleration(
            &dynamics,
            &zero,
            &zero,
            &force,
            &mass_preconditioner,
            moving::initial_config(),
            INITIAL_WORKSPACE,
            cancelled,
        )
        .map_err(|e| e.to_string())?;
        let state = SecondOrderState::new(0.0, &zero, &zero, &initial.value);
        let mut recording = RecordedStructural::new(
            moving::method(n, moving::RHO),
            &dynamics,
            &state,
            StructuralRecordingConfig {
                steps: moving::STEPS,
                adjoint: moving::adjoint_config(n),
                max_workspace_components: 65_536,
            },
        )
        .map_err(|e| e.to_string())?;
        let forward = recording
            .advance(moving::STEPS, moving::STEPS, cancelled)
            .map_err(|e| e.to_string())?;
        if forward.status != StructuralRecordingStatus::ReachedEnd {
            return Err(format!("incomplete plate trial: {:?}", forward.status));
        }
        let adjoint = dynamics
            .effective_preconditioner(moving::STEP, moving::RHO, cancelled)
            .map_err(|e| e.to_string())?;
        let result = recording
            .pullback_samples(
                &self.indices,
                8,
                &Observations {
                    data: self,
                    offset: point[2],
                },
                &adjoint,
                StructuralReplayBudget {
                    checkpoints: 8,
                    forward_steps: 256,
                },
                cancelled,
            )
            .map_err(|e| e.to_string())?;
        let initial = second_order_acceleration_vjp(
            &dynamics,
            &zero,
            &zero,
            &force,
            &result.gradient.initial_a,
            &mass_preconditioner,
            &mass_preconditioner,
            moving::initial_config(),
            INITIAL_WORKSPACE,
            cancelled,
        )
        .map_err(|e| e.to_string())?;
        let mut initial_load = vec![0.0; 6];
        dynamics.forcing_vjp(0.0, &initial.forcing, &mut initial_load)?;
        let gradient: Vec<_> = [0, 4, 5]
            .iter()
            .map(|&i| result.gradient.parameters[i] + initial.parameters[i] + initial_load[i])
            .collect();
        if cancelled() {
            return Err("plate calibration cancelled".into());
        }
        let mut ci = Vec::with_capacity(6);
        let mut ji = vec![0.0; 18];
        for (i, (x, b)) in point.iter().zip(BOUNDS).enumerate() {
            ci.extend([b[0] - x, x - b[1]]);
            ji[6 * i + i] = -1.0;
            ji[(2 * i + 1) * 3 + i] = 1.0;
        }
        Ok(Some(SqpSample {
            f: result.value,
            gradient,
            ce: Vec::new(),
            ci,
            je: Vec::new(),
            ji,
        }))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let data = PlateCalibration::synthetic()?;
    let mut evaluate = |point: &[f64]| data.evaluate(point, &mut || false);
    let mut study = SqpState::try_new(&START, 9, &mut evaluate, None).map_err(|e| e.to_string())?;
    let initial = study.sample().f;
    let result = study
        .try_run(&mut evaluate, 1e-7, 80, 300, None)
        .map_err(|e| e.to_string())?;
    println!(
        "source=synthetic-dense-DKT model=clamped-plate triangles={} free_dofs={} observations={} stop={:?}",
        data.mesh.tris.len(),
        data.plate.free,
        data.indices.len(),
        result.stop
    );
    println!(
        "stiffness_scale={:.9} load_N={:.9} offset_m={:.12e} mass_scale=1",
        study.point()[0],
        study.point()[1],
        study.point()[2] * DISPLACEMENT_UNIT
    );
    println!(
        "objective_before={initial:.12e} objective_after={:.12e} iterations={} evaluations={}",
        study.sample().f,
        study.iterations(),
        study.evaluations()
    );
    println!(
        "scope=fixed-mesh numerical recovery; physical and numerical damping remain distinct; no experimental validation"
    );
    if result.stop != SqpStop::Converged {
        return Err("plate fit stopped before its local KKT tolerance".into());
    }
    Ok(())
}
