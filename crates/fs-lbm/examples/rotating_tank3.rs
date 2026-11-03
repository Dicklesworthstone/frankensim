//! Smooth spin-up of a closed, fully filled tank on the native D3Q19 solver.
//!
//! Run `cargo run -p fs-lbm --example rotating_tank3 -- 12` (8, 12 or 16 cells).
//! Inputs are lattice units: dx=dt=rho0=1, BGK tau=0.8, nu=(tau-0.5)/3.
//! All six walls are stationary in the tank's co-moving coordinates. The
//! prescribed z rotation ramps smoothly; no gravity or other body force acts.
//! Fluid angular momentum uses the native Guo force-corrected relative velocity.
//! Wall impulse and held frame-source impulse close its measured balance.
//!
//! This is a small weakly compressible flow example with explicit low-Mach
//! inputs, not free-surface pouring, thermal transport, SI calibration, a mesh
//! convergence study or an exact global-energy claim. Vessel motion is prescribed.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, FaceBoundary3, RigidFrame3};
use fs_math::det;
use std::error::Error;

const STEPS: usize = 320;
const SPINUP: f64 = 160.0;
const TAU: f64 = 0.8;
const OMEGA: f64 = 0.002;

/// Quintic ramp of angular velocity, with zero acceleration and jerk at each end.
fn motion(time: f64) -> (f64, f64) {
    let s = (time / SPINUP).min(1.0);
    let omega = OMEGA * s * s * s * (10.0 - 15.0 * s + 6.0 * s * s);
    let alpha = 30.0 * OMEGA / SPINUP * s * s * (1.0 - s) * (1.0 - s);
    (omega, alpha)
}

struct Observation {
    mass: f64,
    angular_momentum_z: f64,
    source_torque_z: f64,
    maximum_mach: f64,
}

fn observe(grid: &BoundaryGrid3, origin: [f64; 3]) -> Observation {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Observation {
        mass: 0.0,
        angular_momentum_z: 0.0,
        source_torque_z: 0.0,
        maximum_mach: 0.0,
    };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let rho = grid.density(x, y, z);
                let u = grid.velocity(x, y, z); // includes the Guo half-force
                let force = grid.force_density(x, y, z);
                let rx = x as f64 + 0.5 - origin[0];
                let ry = y as f64 + 0.5 - origin[1];
                result.mass += rho;
                result.angular_momentum_z += rho * (rx * u[1] - ry * u[0]);
                result.source_torque_z += rx * force[1] - ry * force[0];
                result.maximum_mach = result
                    .maximum_mach
                    .max(det::sqrt(3.0 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2])));
            }
        }
    }
    result
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let n: usize = args.next().as_deref().unwrap_or("12").parse()?;
    if ![8, 12, 16].contains(&n) || args.next().is_some() {
        return Err("usage: rotating_tank3 [8|12|16]".into());
    }
    let origin = [0.5 * n as f64; 3]; // native cell centres are x+0.5
    let walls = BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]);
    let mut grid = BoundaryGrid3::new(n, n, n, TAU, [0.0; 3], walls);
    let links = grid.boundary_links(); // fixed topology; retain every wall link
    let gate = CancelGate::new();
    let initial = observe(&grid, origin);
    let mut wall_impulse_z = 0.0;
    let mut frame_impulse_z = 0.0;
    println!(
        "# {n}^3 closed filled cells; lattice dx=dt=rho0=1; tau={TAU}; nu={:.6}",
        grid.viscosity()
    );
    println!("# prescribed rotation: omega_final={OMEGA}; spinup={SPINUP}; centre={origin:?}");
    println!("# wall_Jz is impulse ON the walls; frame_Jz is the held inertial-source impulse");
    println!(
        "step,omega_mid,alpha_mid,mass_drift_rel,fluid_relative_Lz,wall_Jz,frame_Jz,balance_residual_Lz,max_relative_Mach"
    );
    for step in 0..STEPS {
        let (omega, alpha) = motion(step as f64 + 0.5);
        RigidFrame3 {
            origin,
            angular_velocity: [0.0, 0.0, omega],
            angular_acceleration: [0.0, 0.0, alpha],
            origin_acceleration: [0.0; 3],
        }
        .prepare_step(&mut grid, [0.0; 3], &gate)?;
        // Preparation preserves the previous physical state. Only this actual
        // collision/stream step supplies fluid and wall momentum exchange.
        let exchange = grid.step_with_wall_momentum(&links, origin);
        let state = observe(&grid, origin);
        wall_impulse_z += exchange.wall_angular_impulse[2];
        // The force field is held unchanged during the unit-time native step.
        frame_impulse_z += state.source_torque_z;
        if !state.mass.is_finite()
            || !state.angular_momentum_z.is_finite()
            || !state.maximum_mach.is_finite()
            || state.maximum_mach > 0.1
        {
            return Err("flow left the example's finite, relative Mach <= 0.1 domain".into());
        }
        if step == 0 || (step + 1) % 40 == 0 {
            let residual = state.angular_momentum_z - initial.angular_momentum_z + wall_impulse_z
                - frame_impulse_z;
            println!(
                "{},{omega:.8e},{alpha:.8e},{:.8e},{:.8e},{wall_impulse_z:.8e},{frame_impulse_z:.8e},{residual:.8e},{:.8e}",
                step + 1,
                (state.mass - initial.mass) / initial.mass,
                state.angular_momentum_z,
                state.maximum_mach
            );
        }
    }
    Ok(())
}
