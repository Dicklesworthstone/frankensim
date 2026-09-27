//! A bounded Newton/FGMRES solve of a meshed plate with physical Jacobi scaling.
//!
//! `cargo run -p fs-plate --example preconditioned_plate`
//! `cargo run -p fs-plate --example preconditioned_plate -- 12`
//!
//! The optional mesh size is 8 (147 free DOFs, default) or 12 (363 free DOFs).
//!
//! The identity comparison uses exactly the same mesh, state, force and solver
//! controls. Iteration counts describe this experiment, not mesh-independent
//! convergence or a hardware timing benchmark.
#[path = "moving_plate.rs"]
#[allow(dead_code)]
pub(crate) mod moving;

use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateModel, PlateSection, assemble};
use fs_plate_transient::PlateDynamics;
use fs_solver::{NewtonKrylovConfig, NewtonReport};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, TimeSolveError,
    second_order_adjoint::{SecondOrderVjp, trajectory::StructuralTrajectoryModel},
};
use std::error::Error;

/// Fixed physical inputs: stiffness, mass, two Rayleigh coefficients and force.
pub const PARAMETERS: [f64; 5] = [1.1, 0.9, 3.0, 2e-4, 1.7];
/// Time increment in seconds for the bounded larger-mesh comparison.
pub const STEP: f64 = 0.0002;
/// Generalized-alpha high-frequency spectral radius.
pub const RHO: f64 = 0.6;
/// Maximum basis length retained in each forward FGMRES cycle.
pub const RESTART: usize = 16;
/// Maximum forward FGMRES cycles per Newton attempt.
pub const CYCLES: usize = 8;
/// Maximum Newton attempts per physical step.
pub const NEWTON_ATTEMPTS: usize = 8;

/// Assemble a clamped 1m by 0.6m plate with three DKT DOFs per free node.
pub fn build_plate(elements: usize) -> Result<(PlateMesh, PlateModel), Box<dyn Error>> {
    let mesh = PlateMesh::rectangle(1.0, 0.6, elements, elements);
    let section = PlateSection::isotropic(3.0e9, 0.28, 0.008, 600.0)?;
    let plate = assemble(
        &mesh,
        &section,
        &PlateMesh::rectangle_boundary(elements, elements),
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )?;
    Ok((mesh, plate))
}

/// An explicit short-restart policy; no basis grows with the plate dimension.
#[must_use]
pub fn method(n: usize, step: f64, restart: usize, cycles: usize) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        step,
        RHO,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-11,
                linear_restart: restart,
                max_linear_cycles: cycles,
                forcing_minimum: 1e-12,
                forcing_maximum: 1e-8,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: NEWTON_ATTEMPTS,
        },
    )
}

/// Identical physics with the trait's default identity forward preconditioner.
pub struct IdentityProblem<'a, P: ?Sized>(pub &'a P);

impl<P: SecondOrderProblem + ?Sized> SecondOrderProblem for IdentityProblem<'_, P> {
    fn dimension(&self) -> usize {
        self.0.dimension()
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.0.mass_apply(input, output);
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        self.0.damping_apply(input, output);
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        self.0.internal_force(q, output);
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        self.0.tangent_apply(q, direction, output);
    }
}

impl<P: SecondOrderVjp + ?Sized> SecondOrderVjp for IdentityProblem<'_, P> {
    fn parameter_count(&self) -> usize {
        self.0.parameter_count()
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.0.mass_transpose_apply(seed, output)
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.0.damping_transpose_apply(seed, output)
    }
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.0.tangent_transpose_apply(q, seed, output)
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        parameters: &mut [f64],
    ) -> Result<(), String> {
        self.0.residual_parameter_vjp(q, v, a, seed, parameters)
    }
}

impl<P: StructuralTrajectoryModel + ?Sized> StructuralTrajectoryModel for IdentityProblem<'_, P> {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        self.0.forcing(time, output)
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], parameters: &mut [f64]) -> Result<(), String> {
        self.0.forcing_vjp(time, seed, parameters)
    }
}

fn print_report(label: &str, report: &NewtonReport) {
    let iterations: usize = report
        .history
        .iter()
        .map(|step| step.linear_iterations)
        .sum();
    if report.converged {
        println!(
            "{label}: converged=true; inner_iterations={iterations}; true_newton_residual={:.12e}",
            report.residual_norm
        );
    } else {
        // A failed inner solve is not retained in Newton's accepted-attempt
        // history. Do not present that partial history as its iteration count.
        println!(
            "{label}: converged=false; true_newton_residual={:.12e}; diagnosis={:?}",
            report.residual_norm, report.diagnosis
        );
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = std::env::args().skip(1);
    let elements = match arguments.next().as_deref() {
        None | Some("8") => 8,
        Some("12") => 12,
        Some(_) => return Err("expected mesh size 8 or 12".into()),
    };
    if arguments.next().is_some() {
        return Err("usage: preconditioned_plate [8|12]".into());
    }
    let (mesh, plate) = build_plate(elements)?;
    let load = moving::MovingLoad {
        mesh: &mesh,
        plate: &plate,
        amplitude: PARAMETERS[4],
    };
    let dynamics = PlateDynamics::new(
        &plate,
        moving::dynamics_parameters(PARAMETERS),
        &load,
        moving::dynamics_budget(&plate),
        &mut || false,
    )?;
    let initial = moving::resting_state(&dynamics)?;
    let mut preconditioned = initial.clone();
    let solver = method(plate.free, STEP, RESTART, CYCLES);
    let mut force = vec![0.0; plate.free];
    dynamics.forcing(solver.forcing_time(initial.t)?, &mut force)?;
    println!(
        "clamped_plate: triangles={}; dofs={}; step_s={STEP}; restart={RESTART}; cycles_per_newton={CYCLES}; newton_attempt_cap={NEWTON_ATTEMPTS}; inner_iteration_cap_per_newton={}",
        mesh.tris.len(),
        plate.free,
        RESTART * CYCLES
    );
    let report = solver.step(&mut preconditioned, &dynamics, &force)?;
    print_report("physical_jacobi", &report.newton);
    let mut baseline = initial;
    match solver.step(&mut baseline, &IdentityProblem(&dynamics), &force) {
        Ok(report) => print_report("identity", &report.newton),
        Err(TimeSolveError::NotConverged(report)) => print_report("identity", &report),
        Err(error) => return Err(error.into()),
    }
    let weights = moving::sensor_weights(&mesh, &plate)?;
    let displacement: f64 = weights
        .iter()
        .zip(&preconditioned.q)
        .map(|(w, q)| w * q)
        .sum();
    println!("endpoint_sensor_deflection_m={displacement:.12e}");
    Ok(())
}
