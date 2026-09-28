//! A moving transverse force on a meshed clamped plate, with sampled adjoints.
//!
//! `cargo run -p fs-plate-transient --example moving_plate`
//!
//! Geometry and material data are explicit illustrative SI inputs. This uses
//! all assembled DKT displacement/slope DOFs and the shared time integrator.
//! It demonstrates a discrete design derivative, not experimental validation.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateMesh, PlateModel, PlateSection, assemble,
    loading::{PlateLoadBudget, PlatePointStencil},
};
use fs_plate_transient::{PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateLoad};
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
            RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
            StructuralReplayBudget, StructuralTrajectoryModel, samples::StructuralSampleObjective,
        },
    },
};
use std::error::Error;

/// Fixed step in seconds; observation indices never interpolate the solution.
pub const STEP: f64 = 0.002;
/// Length of the forward trajectory.
pub const STEPS: usize = 12;
/// High-frequency spectral radius for generalized-alpha.
pub const RHO: f64 = 0.6;
/// Initial and three subsequent sensor observations.
pub const SAMPLES: [usize; 4] = [0, 4, 8, 12];
/// Displacement scale used only to nondimensionalize the objective.
pub const DEFLECTION_SCALE_M: f64 = 1.0e-4;

/// A 1m by 0.6m clamped plate: 32 DKT triangles and 27 free DOFs.
pub fn build_plate() -> Result<(PlateMesh, PlateModel), Box<dyn Error>> {
    let mesh = PlateMesh::rectangle(1.0, 0.6, 4, 4);
    let section = PlateSection::isotropic(3.0e9, 0.28, 0.008, 600.0)?;
    let plate = assemble(
        &mesh,
        &section,
        &PlateMesh::rectangle_boundary(4, 4),
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )?;
    Ok((mesh, plate))
}

/// Mechanical entries precede the moving force amplitude, in the adapter order.
#[must_use]
pub fn dynamics_parameters(p: [f64; 5]) -> PlateDynamicsParameters {
    PlateDynamicsParameters {
        stiffness_scale: p[0],
        mass_scale: p[1],
        mass_damping_per_s: p[2],
        stiffness_damping_s: p[3],
    }
}
/// Admit exactly this spatial model and one applied-load parameter.
#[must_use]
pub fn dynamics_budget(plate: &PlateModel) -> PlateDynamicsBudget {
    PlateDynamicsBudget {
        max_dofs: plate.free,
        max_nonzeros: plate.k.nnz() + plate.m.nnz(),
        max_full_dofs: plate.dof_map.len(),
        max_load_parameters: 1,
    }
}

/// Constant-amplitude force travelling across the fixed plate, in newtons.
/// The path is prescribed; its coordinates are not design parameters.
pub struct MovingLoad<'a> {
    /// Fixed topology used to locate the physical force point.
    pub mesh: &'a PlateMesh,
    /// Assembled support map and reduced force-vector dimension.
    pub plate: &'a PlateModel,
    /// Transverse force amplitude in newtons.
    pub amplitude: f64,
}
impl MovingLoad<'_> {
    /// Locate the actual force point without nearest-node snapping.
    pub fn stencil(&self, time: f64) -> Result<PlatePointStencil, String> {
        if !time.is_finite() || !(0.0..=STEP * STEPS as f64 + 1e-12).contains(&time) {
            return Err("moving-load time is outside its prescribed path".into());
        }
        let point = [0.12 + 0.76 * time / (STEP * STEPS as f64), 0.27];
        PlatePointStencil::locate(
            self.mesh,
            self.plate,
            point,
            PlateLoadBudget {
                max_nodes: self.mesh.nodes.len(),
                max_triangles: self.mesh.tris.len(),
            },
        )
        .map_err(|e| e.to_string())
    }
}
impl PlateLoad for MovingLoad<'_> {
    fn parameter_count(&self) -> usize {
        1
    }
    fn forcing(&self, time: f64, out: &mut [f64]) -> Result<(), String> {
        out.fill(0.0);
        self.stencil(time)?
            .add_load([self.amplitude, 0.0, 0.0], out)
            .map_err(|e| e.to_string())?;
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], parameters: &mut [f64]) -> Result<(), String> {
        if parameters.len() != 1 {
            return Err("one moving-force amplitude parameter required".into());
        }
        parameters[0] = self
            .stencil(time)?
            .load_vjp(seed)
            .map_err(|e| e.to_string())?[0];
        Ok(())
    }
}

/// P1 displacement sensor at a fixed physical point inside the plate.
pub fn sensor_weights(mesh: &PlateMesh, plate: &PlateModel) -> Result<Vec<f64>, Box<dyn Error>> {
    let sensor = PlatePointStencil::locate(
        mesh,
        plate,
        [0.55, 0.32],
        PlateLoadBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        },
    )?;
    let mut weights = vec![0.0; plate.free];
    sensor.add_load([1.0, 0.0, 0.0], &mut weights)?;
    Ok(weights)
}

/// Production Newton/Krylov configuration, retaining every spatial DOF.
#[must_use]
pub fn method(n: usize, rho: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        STEP,
        rho,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-11,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_minimum: 1e-12,
                forcing_maximum: 1e-10,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 8,
        },
    )
}
/// Explicit transposed-solve controls; Jacobi is supplied separately.
#[must_use]
pub fn adjoint_config(n: usize) -> SecondOrderAdjointConfig {
    SecondOrderAdjointConfig {
        restart: n,
        max_cycles: 4,
        tolerance: 1e-12,
    }
}
/// Diagonal mass solves converge with the admitted exact mass preconditioner.
#[must_use]
pub fn initial_config() -> InitialSolveConfig {
    InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1e-13,
    }
}

/// Start at zero displacement/velocity with a dynamically consistent acceleration.
pub fn resting_state<L: PlateLoad + ?Sized>(
    dynamics: &PlateDynamics<'_, L>,
) -> Result<SecondOrderState, Box<dyn Error>> {
    let n = dynamics.model().free;
    let zero = vec![0.0; n];
    let mut force = vec![0.0; n];
    dynamics.forcing(0.0, &mut force)?;
    let preconditioner = dynamics.mass_preconditioner();
    let acceleration = second_order_acceleration(
        dynamics,
        &zero,
        &zero,
        &force,
        &preconditioner,
        initial_config(),
        initial_workspace_components(n, dynamics.parameter_count(), initial_config())
            .ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    Ok(SecondOrderState::new(
        0.0,
        &zero,
        &zero,
        &acceleration.value,
    ))
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
struct DeflectionObjective<'a>(&'a [f64]);
impl StructuralSampleObjective for DeflectionObjective<'_> {
    fn evaluate(
        &self,
        _sample: usize,
        state: &SecondOrderState,
        bars: (&mut [f64], &mut [f64], &mut [f64]),
        parameters: &mut [f64],
    ) -> Result<f64, String> {
        let value = dot(self.0, &state.q) / DEFLECTION_SCALE_M;
        for (bar, weight) in bars.0.iter_mut().zip(self.0) {
            *bar = value * weight / DEFLECTION_SCALE_M;
        }
        bars.1.fill(0.0);
        bars.2.fill(0.0);
        parameters.fill(0.0);
        Ok(0.5 * value * value)
    }
}

/// Measured endpoint displacements and the total physical-parameter gradient.
pub struct PlateRun {
    /// Accepted observation times in seconds and sensor displacement in metres.
    pub samples: Vec<(f64, f64)>,
    /// Sum of squared normalized sensor displacements, divided by two.
    pub loss: f64,
    /// Derivatives in the four dynamics parameters, then force amplitude.
    pub gradient: Vec<f64>,
    /// All free displacement, slope, velocity and acceleration components.
    pub final_state: SecondOrderState,
    /// Complete forward steps used by the checkpointed reverse sweep.
    pub replays: usize,
    /// Maximum simultaneously parked q/v/a checkpoints.
    pub checkpoints: usize,
}

/// One complete forward/adjoint experiment, including initial force dependence.
pub fn evaluate(
    mesh: &PlateMesh,
    plate: &PlateModel,
    p: [f64; 5],
) -> Result<PlateRun, Box<dyn Error>> {
    let load = MovingLoad {
        mesh,
        plate,
        amplitude: p[4],
    };
    let dynamics = PlateDynamics::new(
        plate,
        dynamics_parameters(p),
        &load,
        dynamics_budget(plate),
        &mut || false,
    )?;
    let initial = resting_state(&dynamics)?;
    let method = method(plate.free, RHO);
    let adjoint = adjoint_config(plate.free);
    let mut tape = RecordedStructural::new(
        method,
        &dynamics,
        &initial,
        StructuralRecordingConfig {
            steps: STEPS,
            adjoint,
            max_workspace_components: method
                .adjoint_workspace_components(5, adjoint)
                .ok_or("trajectory workspace overflow")?,
        },
    )?;
    let weights = sensor_weights(mesh, plate)?;
    let mut samples = Vec::with_capacity(SAMPLES.len());
    for endpoint in SAMPLES {
        let report = tape.advance(endpoint - tape.accepted_steps(), STEPS, &mut || false)?;
        if endpoint == STEPS && report.status != StructuralRecordingStatus::ReachedEnd {
            return Err("moving-plate forward trajectory incomplete".into());
        }
        samples.push((tape.time(), dot(&weights, &tape.state().q)));
    }
    let jacobi = dynamics.effective_preconditioner(STEP, RHO, &mut || false)?;
    let result = tape.pullback_samples(
        &SAMPLES,
        SAMPLES.len(),
        &DeflectionObjective(&weights),
        &jacobi,
        StructuralReplayBudget {
            checkpoints: tape.required_checkpoints(),
            forward_steps: 128,
        },
        &mut || false,
    )?;
    let mass = dynamics.mass_preconditioner();
    let mut force = vec![0.0; plate.free];
    dynamics.forcing(0.0, &mut force)?;
    let initial_gradient = second_order_acceleration_vjp(
        &dynamics,
        &initial.q,
        &initial.v,
        &force,
        &result.gradient.initial_a,
        &mass,
        &mass,
        initial_config(),
        initial_workspace_components(plate.free, 5, initial_config())
            .ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    let mut force_gradient = vec![0.0; 5];
    dynamics.forcing_vjp(0.0, &initial_gradient.forcing, &mut force_gradient)?;
    let mut gradient = result.gradient.parameters;
    for i in 0..5 {
        gradient[i] += initial_gradient.parameters[i] + force_gradient[i];
    }
    Ok(PlateRun {
        samples,
        loss: result.value,
        gradient,
        final_state: tape.state().clone(),
        replays: result.gradient.replayed_steps,
        checkpoints: result.gradient.peak_checkpoints,
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let (mesh, plate) = build_plate()?;
    let p = [1.1, 0.9, 3.0, 2e-4, 1.7];
    let run = evaluate(&mesh, &plate, p)?;
    println!(
        "clamped plate: {} triangles, {} free DKT degrees of freedom",
        mesh.tris.len(),
        plate.free
    );
    println!("time_s,sensor_deflection_m");
    for (time, deflection) in &run.samples {
        println!("{time:.6},{deflection:.12e}");
    }
    println!(
        "sampled_loss={:.12e}; checkpoints={}; forward_replays={}",
        run.loss, run.checkpoints, run.replays
    );
    for (name, gradient) in [
        "stiffness_scale",
        "mass_scale",
        "mass_damping_per_s",
        "stiffness_damping_s",
        "force_amplitude_n",
    ]
    .iter()
    .zip(&run.gradient)
    {
        println!("d_loss/d_{name}={gradient:.12e}");
    }
    // One prescribed trial verifies a descent direction; this is not a fit or
    // an optimizer. Common scaling of K/M/F is unidentifiable in this experiment.
    let scales = [1.0, 1.0, 10.0, 1e-3, 1.0];
    let norm = run
        .gradient
        .iter()
        .zip(scales)
        .map(|(g, s)| (g * s).powi(2))
        .sum::<f64>()
        .sqrt();
    if !(norm.is_finite() && norm > 0.0) {
        return Err("no finite descent direction".into());
    }
    let candidate =
        std::array::from_fn(|i| p[i] - 0.001 * scales[i] * scales[i] * run.gradient[i] / norm);
    let trial = evaluate(&mesh, &plate, candidate)?;
    if trial.loss >= run.loss {
        return Err("prescribed gradient trial did not reduce sampled loss".into());
    }
    println!(
        "one_scaled_gradient_trial_loss={:.12e}; decrease={:.6}%",
        trial.loss,
        100.0 * (run.loss - trial.loss) / run.loss
    );
    Ok(())
}
