//! Fit one authored Hertz stiffness to synthetic velocity observations on a
//! spatial DKT plate: `cargo run -p fs-plate-transient --example calibrate_plate_impact`.
//!
//! This is a bounded local inverse problem with all other controls fixed. Its
//! observations come from the same numerical model, not experimental evidence.
//! Material-source contacts require re-resolution and cannot opt into these
//! independently authored contact parameters. All gradients below include the
//! consistent initial acceleration and the original parameterized plate load.

#[path = "moving_plate.rs"]
#[allow(dead_code)]
mod spatial_plate;

use fs_dcontact::Obstacle;
use fs_plate::{
    PlateMesh, PlateModel,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact, PlateImpactBudget,
    PlateLoad,
};
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState,
    initialization::{
        InitialSolveConfig, initial_workspace_components, second_order_acceleration,
        second_order_acceleration_vjp,
    },
    second_order_adjoint::{
        SecondOrderAdjointConfig, SecondOrderVjp,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralReplayBudget,
            StructuralTrajectoryModel, samples::StructuralSampleObjective,
        },
    },
};
use std::error::Error;

/// Plate scales/damping, load amplitude, then body mass, contact K, chi and gap.
pub type Parameters = [f64; 9];
/// Observation endpoints on the fixed eight-step trajectory.
pub const SAMPLES: [usize; 3] = [2, 4, 8];
const VELOCITY_SCALE: f64 = 0.1;

/// A real independently parameterized applied load preceding contact parameters.
pub struct DrivingLoad {
    /// Applied force amplitude in newtons.
    pub amplitude: f64,
    weights: Vec<f64>,
}
impl PlateLoad for DrivingLoad {
    fn parameter_count(&self) -> usize {
        1
    }
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        if output.len() != self.weights.len() || !time.is_finite() {
            return Err("invalid calibration driving load".into());
        }
        let factor = self.amplitude * (1.0 + 20.0 * time);
        for (out, weight) in output.iter_mut().zip(&self.weights) {
            *out = factor * weight;
        }
        if output.iter().all(|x| x.is_finite()) {
            Ok(())
        } else {
            Err("nonfinite calibration driving load".into())
        }
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], parameters: &mut [f64]) -> Result<(), String> {
        if seed.len() != self.weights.len() || parameters.len() != 1 || !time.is_finite() {
            return Err("invalid calibration load derivative".into());
        }
        parameters[0] = (1.0 + 20.0 * time)
            * self
                .weights
                .iter()
                .zip(seed)
                .map(|(a, b)| a * b)
                .sum::<f64>();
        if parameters[0].is_finite() {
            Ok(())
        } else {
            Err("nonfinite load derivative".into())
        }
    }
}

fn patch(mesh: &PlateMesh, model: &PlateModel) -> Result<PlatePatchStencil, Box<dyn Error>> {
    Ok(PlatePatchStencil::locate(
        mesh,
        model,
        &[
            PlatePatchPoint {
                point: [0.4, 0.2],
                area_m2: 0.001,
            },
            PlatePatchPoint {
                point: [0.6, 0.2],
                area_m2: 0.003,
            },
        ],
        PlatePatchBudget {
            max_points: 2,
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    )?)
}

/// Use the same footprint for the independently prescribed plate force.
pub fn driving_load(
    mesh: &PlateMesh,
    plate: &PlateModel,
    amplitude: f64,
) -> Result<DrivingLoad, Box<dyn Error>> {
    let mut weights = vec![0.0; plate.free];
    patch(mesh, plate)?.add_load([1.0, 0.0, 0.0], &mut weights)?;
    Ok(DrivingLoad { amplitude, weights })
}

/// Build the authored contact without opting into calibration yet.
pub fn authored_impact<'a>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    p: Parameters,
    load: &'a DrivingLoad,
) -> Result<PlateImpact<'a, DrivingLoad>, Box<dyn Error>> {
    let dynamics = PlateDynamics::new(
        plate,
        PlateDynamicsParameters {
            stiffness_scale: p[0],
            mass_scale: p[1],
            mass_damping_per_s: p[2],
            stiffness_damping_s: p[3],
        },
        load,
        PlateDynamicsBudget {
            max_dofs: plate.free,
            max_nonzeros: plate.k.nnz() + plate.m.nnz(),
            max_full_dofs: plate.dof_map.len(),
            max_load_parameters: 1,
        },
        &mut || false,
    )?;
    let law = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![p[8]],
        vec![1.0],
        p[6],
        1.5,
        "authored synthetic plate-impact calibration".into(),
    )?
    .with_internal_loss(p[7])?;
    Ok(PlateImpact::new(
        dynamics,
        patch(mesh, plate)?,
        p[5],
        law,
        PlateImpactBudget {
            max_free_dofs: plate.free + 1,
            max_contact_points: 1,
        },
    )?)
}

/// Shared bounded Newton solve, with fixed time discretization during fitting.
pub fn method(n: usize, step: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        step,
        0.6,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1.0e-12,
                relative_tolerance: 1.0e-12,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_minimum: 1.0e-13,
                forcing_maximum: 1.0e-10,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 16,
        },
    )
}
/// Explicit bounded transpose solve.
pub fn adjoint_config(n: usize) -> SecondOrderAdjointConfig {
    SecondOrderAdjointConfig {
        restart: n,
        max_cycles: 4,
        tolerance: 1.0e-12,
    }
}
fn initial_config() -> InitialSolveConfig {
    InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1.0e-12,
    }
}

/// Recompute the physical initial acceleration at every parameter point.
pub fn consistent_state(
    model: &PlateImpact<'_, DrivingLoad>,
    q: &[f64],
    v: &[f64],
) -> Result<SecondOrderState, Box<dyn Error>> {
    let mut force = vec![0.0; model.dimension()];
    model.forcing(0.0, &mut force)?;
    let config = initial_config();
    let initial = second_order_acceleration(
        model,
        q,
        v,
        &force,
        &model.mass_preconditioner(),
        config,
        initial_workspace_components(model.dimension(), model.parameter_count(), config)
            .ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    Ok(SecondOrderState::new(0.0, q, v, &initial.value))
}

/// Chain the initial-acceleration seed through contact, mass and applied forcing.
pub fn add_initial_pullback(
    model: &PlateImpact<'_, DrivingLoad>,
    initial: &SecondOrderState,
    acceleration_seed: &[f64],
    parameters: &mut [f64],
) -> Result<(), Box<dyn Error>> {
    let mut force = vec![0.0; model.dimension()];
    model.forcing(0.0, &mut force)?;
    let mass = model.mass_preconditioner();
    let config = initial_config();
    let result = second_order_acceleration_vjp(
        model,
        &initial.q,
        &initial.v,
        &force,
        acceleration_seed,
        &mass,
        &mass,
        config,
        initial_workspace_components(model.dimension(), model.parameter_count(), config)
            .ok_or("initial VJP workspace overflow")?,
        &mut || false,
    )?;
    let mut applied = vec![0.0; model.parameter_count()];
    model.forcing_vjp(0.0, &result.forcing, &mut applied)?;
    for ((out, initial), load) in parameters.iter_mut().zip(result.parameters).zip(applied) {
        *out += initial + load;
    }
    Ok(())
}

struct VelocityFit<'a> {
    body: usize,
    targets: &'a [f64; 3],
}
impl StructuralSampleObjective for VelocityFit<'_> {
    fn evaluate(
        &self,
        sample: usize,
        state: &SecondOrderState,
        bars: (&mut [f64], &mut [f64], &mut [f64]),
        parameters: &mut [f64],
    ) -> Result<f64, String> {
        let residual = (state.v[self.body] - self.targets[sample]) / VELOCITY_SCALE;
        bars.0.fill(0.0);
        bars.1.fill(0.0);
        bars.2.fill(0.0);
        parameters.fill(0.0);
        bars.1[self.body] = residual / VELOCITY_SCALE;
        Ok(0.5 * residual * residual)
    }
}
/// The actual numerical observations and optional complete parameter gradient.
pub struct FitRun {
    /// Body velocity observations in metres per second.
    pub samples: [f64; 3],
    /// Sum of squared normalized observation residuals divided by two.
    pub loss: f64,
    /// Original five plate/load entries then mass, K, chi and gap; absent on forward trials.
    pub gradient: Option<Vec<f64>>,
}

/// Record/replay the original shared trajectory; no solver iterations are differentiated.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(
    mesh: &PlateMesh,
    plate: &PlateModel,
    p: Parameters,
    q: &[f64],
    v: &[f64],
    step: f64,
    targets: &[f64; 3],
    differentiate: bool,
) -> Result<FitRun, Box<dyn Error>> {
    let load = driving_load(mesh, plate, p[4])?;
    let model = authored_impact(mesh, plate, p, &load)?.with_contact_calibration()?;
    let initial = consistent_state(&model, q, v)?;
    let method = method(model.dimension(), step);
    let adjoint = adjoint_config(model.dimension());
    let mut tape = RecordedStructural::new(
        method,
        &model,
        &initial,
        StructuralRecordingConfig {
            steps: 8,
            adjoint,
            max_workspace_components: method
                .adjoint_workspace_components(model.parameter_count(), adjoint)
                .ok_or("trajectory workspace overflow")?,
        },
    )?;
    let mut samples = [0.0; 3];
    for (i, endpoint) in SAMPLES.iter().copied().enumerate() {
        tape.advance(endpoint - tape.accepted_steps(), 8, &mut || false)?;
        samples[i] = tape.state().v[plate.free];
    }
    let loss = samples
        .iter()
        .zip(targets)
        .map(|(a, b)| 0.5 * ((a - b) / VELOCITY_SCALE).powi(2))
        .sum();
    if !differentiate {
        return Ok(FitRun {
            samples,
            loss,
            gradient: None,
        });
    }
    let result = tape.pullback_samples(
        &SAMPLES,
        SAMPLES.len(),
        &VelocityFit {
            body: plate.free,
            targets,
        },
        &model.mass_preconditioner(),
        StructuralReplayBudget {
            checkpoints: tape.required_checkpoints(),
            forward_steps: 64,
        },
        &mut || false,
    )?;
    let mut gradient = result.gradient.parameters;
    add_initial_pullback(&model, &initial, &result.gradient.initial_a, &mut gradient)?;
    Ok(FitRun {
        samples,
        loss: result.value,
        gradient: Some(gradient),
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let (mesh, plate) = spatial_plate::build_plate()?;
    let truth: Parameters = [1.0, 1.0, 0.1, 1.0e-5, 0.01, 0.02, 2.0e6, 1.0, 5.0e-5];
    let mut q = vec![0.0; plate.free + 1];
    let mut v = q.clone();
    q[plate.free] = 1.5e-4;
    v[plate.free] = 0.05;
    let step = 5.0e-5;
    let targets = evaluate(&mesh, &plate, truth, &q, &v, step, &[0.0; 3], false)?.samples;
    let mut p = truth;
    p[6] *= 0.7;
    println!("iteration,contact_k_n_per_m_3_2,sampled_loss");
    for iteration in 0..14 {
        let current = evaluate(&mesh, &plate, p, &q, &v, step, &targets, true)?;
        println!("{iteration},{:.12e},{:.12e}", p[6], current.loss);
        if current.loss < 1.0e-12 {
            break;
        }
        let log_gradient = p[6] * current.gradient.as_ref().ok_or("missing gradient")?[6];
        if !log_gradient.is_finite() || log_gradient == 0.0 {
            break;
        }
        let mut accepted = false;
        let mut fraction = 0.25;
        for _ in 0..12 {
            let mut candidate = p;
            candidate[6] *= (-fraction * log_gradient.signum()).exp();
            let trial = evaluate(&mesh, &plate, candidate, &q, &v, step, &targets, false)?;
            if trial.loss <= current.loss - 1.0e-4 * fraction * log_gradient.abs() {
                p = candidate;
                accepted = true;
                break;
            }
            fraction *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    let final_run = evaluate(&mesh, &plate, p, &q, &v, step, &targets, false)?;
    println!(
        "target_K={:.12e}; fitted_K={:.12e}; final_loss={:.12e}",
        truth[6], p[6], final_run.loss
    );
    Ok(())
}
