//! Native spin-down of a finite-inertia disk against a closed viscous fluid.
//!
//! Run `cargo run -p fs-lbm --example finite_rotor3 -- 12` (8, 12 or 16 cells).
//! The lattice and housing stay in an inertial frame. Fluid traction changes
//! the disk's actual angular velocity at every step; there is no prescribed
//! spin schedule, external motor torque, gravity or added drag coefficient.
//!
//! Inputs use lattice units dx=dt=rho0=1. The disk is an admitted patch of the
//! bottom wall, with ideal bearings constraining its centre and spin axis.
//! Its inertia is that of a uniform thin disk with declared areal mass. The
//! native discretized patch approximates the circle; this is an illustrative
//! low-Mach calculation, not an SI-calibrated or mesh-converged device.
//! Rotor work is reported separately from viscous/numerical fluid dissipation.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3, WallRotor3};
use fs_math::det;
use std::error::Error;

const STEPS: usize = 320;
const TAU: f64 = 0.8;
const INITIAL_SPIN: f64 = 0.003;
const DISK_AREAL_MASS: f64 = 5.0;

struct Observation {
    mass: f64,
    angular_momentum_z: f64,
    maximum_mach: f64,
}

fn observe(grid: &BoundaryGrid3, origin: [f64; 3]) -> Observation {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Observation {
        mass: 0.0,
        angular_momentum_z: 0.0,
        maximum_mach: 0.0,
    };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let rho = grid.density(x, y, z);
                let u = grid.velocity(x, y, z);
                let rx = x as f64 + 0.5 - origin[0];
                let ry = y as f64 + 0.5 - origin[1];
                result.mass += rho;
                result.angular_momentum_z += rho * (rx * u[1] - ry * u[0]);
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
        return Err("usage: finite_rotor3 [8|12|16]".into());
    }
    let center = [0.5 * n as f64, 0.5 * n as f64, 0.0];
    let radius = 0.28 * n as f64;
    let inertia = 0.5 * std::f64::consts::PI * DISK_AREAL_MASS * radius.powi(4);
    let walls = BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]);
    let mut grid = BoundaryGrid3::new(n, n, n, TAU, [0.0; 3], walls);
    let gate = CancelGate::new();
    let mut rotor = WallRotor3::new(
        &grid,
        Face3::ZMin,
        center,
        radius,
        inertia,
        INITIAL_SPIN,
        n * n,
        &gate,
    )?;
    let initial = observe(&grid, center);
    let mut housing_impulse_z = 0.0;
    let mut bearing_impulse_z = 0.0;
    let mut source_impulse_z = 0.0;
    let mut rotor_work_residual = 0.0;
    println!(
        "# {n}^3 closed cells; lattice dx=dt=rho0=1; tau={TAU}; nu={:.6}",
        grid.viscosity()
    );
    println!(
        "# finite disk: radius={radius}; inertia={inertia}; initial_spin={INITIAL_SPIN}; patch_cells={}",
        rotor.patch_cell_count()
    );
    println!(
        "# housing impulses act on the housing; bearing impulses act on the rotor; origin=disk centre"
    );
    println!("# rotor kinetic energy is not the total fluid-plus-rotor energy");
    println!(
        "step,rotor_spin,rotor_angle,rotor_kinetic_energy,mass_drift_rel,fluid_Lz,housing_Jz,bearing_Jz,source_Jz,balance_residual_Lz,rotor_work_residual,max_Mach"
    );
    for step in 0..STEPS {
        let exchange = grid.step_with_finite_inertia_rotor(&mut rotor, 0.0, center, &gate)?;
        housing_impulse_z += exchange.stationary_walls.wall_angular_impulse[2];
        bearing_impulse_z += exchange.bearing_angular_impulse[2];
        source_impulse_z +=
            exchange.external_angular_impulse + exchange.body_force_angular_impulse[2];
        rotor_work_residual += exchange.mechanical_work_residual;
        let state = observe(&grid, center);
        let spin = rotor.angular_velocity();
        let residual = state.angular_momentum_z - initial.angular_momentum_z
            + inertia * (spin - INITIAL_SPIN)
            + housing_impulse_z
            - bearing_impulse_z
            - source_impulse_z;
        if !state.mass.is_finite()
            || !residual.is_finite()
            || !spin.is_finite()
            || !state.maximum_mach.is_finite()
            || state.maximum_mach > 0.1
        {
            return Err("flow left the example's finite, Mach <= 0.1 domain".into());
        }
        if step == 0 || (step + 1) % 40 == 0 {
            println!(
                "{},{spin:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{housing_impulse_z:.8e},{bearing_impulse_z:.8e},{source_impulse_z:.8e},{residual:.8e},{rotor_work_residual:.8e},{:.8e}",
                step + 1,
                rotor.angle(),
                0.5 * inertia * spin * spin,
                (state.mass - initial.mass) / initial.mass,
                state.angular_momentum_z,
                state.maximum_mach
            );
        }
    }
    Ok(())
}
