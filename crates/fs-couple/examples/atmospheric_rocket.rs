//! Source-bearing gas, off-design nozzle, native attitude and atmospheric flight.
//!
//! Run `cargo run -p fs-couple --example atmospheric_rocket --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv 1000 1000000`.
//! See atmospheric_rocket/README.md for the explicit demonstration assumptions.

#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use fs_couple::rocket_blowdown::{
    environment::spherical::SphericalStandardAtmosphere,
    flight::{
        RocketFlightForcing, RocketFlightGeometry, RocketFlightPropulsion, RocketFlightStep,
        advance_atmospheric_rocket,
    },
    rigid::{GasRigidRocketSolve, GasRigidRocketState},
};
use fs_material::{
    gas::atmosphere::{USSA_EFFECTIVE_RADIUS_M, USSA_SURFACE_GRAVITY_M_S2},
    thermochemical_gas::energy::{
        GasEnergySolve,
        entropy::GasEntropySolve,
        parcel::{
            ClosedGasParcel,
            transfer::nozzle::{
                blowdown::{GasBlowdownSolve, RocketBlowdownSolve},
                rocket::{RocketNozzle, RocketNozzleSolve},
            },
        },
    },
};
use fs_mbd::{
    UnitQuaternion, Vec3,
    rocket::{
        RocketState,
        rigid::{RigidRocketInertia, RigidRocketSolve},
    },
};
use std::{
    error::Error,
    io::{Read, Write},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const USAGE: &str = "usage: atmospheric_rocket SOURCE.tsv TEMPERATURE_K PRESSURE_PA [steps=200] [dt_s=0.01] [powered_steps=100]";
const FRAME: u64 = 1;
const DRY_MASS_KG: f64 = 0.1;
const DRAG_COEFFICIENT: f64 = 0.5;
const TANK_VOLUME_M3: f64 = 0.002;
const TANK_LENGTH_M: f64 = 0.2;

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, increment: f64) {
        let corrected = increment - self.correction;
        let next = self.value + corrected;
        self.correction = (next - self.value) - corrected;
        self.value = next;
    }
}

struct Journey {
    state: GasRigidRocketState,
    axial_impulse_n_s: f64,
    expelled_mass_kg: f64,
    mass_residual_kg: f64,
    maximum_mass_residual_kg: f64,
    maximum_axial_impulse_residual_n_s: f64,
    maximum_angular_momentum_residual_kg_m2_s: f64,
    coast_start_mass_kg: Option<f64>,
}

fn controls() -> GasRigidRocketSolve {
    GasRigidRocketSolve {
        gas: RocketBlowdownSolve {
            nozzle: RocketNozzleSolve {
                nozzle: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                relative_mass_flow_tolerance: 1e-9,
            },
            reservoir: GasBlowdownSolve {
                entropy: GasEntropySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    entropy_tolerance_j_kg_k: 1e-10,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-13,
            },
            maximum_donor_fraction: 0.05,
            maximum_iterations: 32,
        },
        mechanics: RigidRocketSolve {
            maximum_substep_s: 0.005,
            maximum_substeps: 10000,
        },
    }
}

fn atmosphere() -> Result<SphericalStandardAtmosphere> {
    // An explicitly declared spherical demonstration, with a 24-hour spin.
    // mu is calibrated to the profile's reference surface g, not a measured
    // Earth gravitational parameter or an ellipsoidal geopotential model.
    let radius = USSA_EFFECTIVE_RADIUS_M;
    Ok(SphericalStandardAtmosphere::new(
        FRAME,
        Vec3::ZERO,
        radius,
        USSA_SURFACE_GRAVITY_M_S2 * radius * radius,
        Vec3::new(0.0, 0.0, core::f64::consts::TAU / 86400.0),
        Vec3::new(0.0, 5.0, 0.0),
    )?)
}

fn initial(reservoir: ClosedGasParcel) -> Result<GasRigidRocketState> {
    let position = Vec3::new(USSA_EFFECTIVE_RADIUS_M + 100.0, 0.0, 0.0);
    let spin = Vec3::new(0.0, 0.0, core::f64::consts::TAU / 86400.0);
    let velocity = spin.cross(position).add(Vec3::new(30.0, 0.0, 0.0));
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.02)?;
    // Uniform gas in the same 2 L, 0.2 m cylindrical tank. The reference
    // mass-distribution model keeps its common COM/principal axes fixed.
    let tank_radius_squared = TANK_VOLUME_M3 / (core::f64::consts::PI * TANK_LENGTH_M);
    let transverse_specific_inertia =
        (3.0 * tank_radius_squared + TANK_LENGTH_M * TANK_LENGTH_M) / 12.0;
    let inertia = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(0.0002, 0.001, 0.001),
        propellant_specific_principal_m2: Vec3::new(
            0.5 * tank_radius_squared,
            transverse_specific_inertia,
            transverse_specific_inertia,
        ),
    };
    let properties = inertia.mass_properties(RocketState {
        position_m: position,
        velocity_m_s: velocity,
        dry_mass_kg: DRY_MASS_KG,
        propellant_mass_kg: reservoir.mass_kg(),
    })?;
    let angular_momentum = Vec3::new(0.0, 0.0, 0.02 * properties.principal_inertia_body().z);
    Ok(GasRigidRocketState::new(
        reservoir,
        DRY_MASS_KG,
        position,
        velocity,
        orientation,
        angular_momentum,
        inertia,
        FRAME,
    )?)
}

fn geometry() -> RocketFlightGeometry {
    RocketFlightGeometry {
        aerodynamic_application_arm_body_m: Vec3::new(-0.04, 0.0, 0.0),
        aerodynamic_reference_length_m: 0.2,
        aerodynamic_reference_area_m2: 0.003,
        nozzle_ambient_arm_body_m: Vec3::new(-0.1, 0.0, 0.0),
    }
}

fn forcing(powered: bool) -> RocketFlightForcing {
    let exit_area = 2e-6;
    let exit_radius_squared = exit_area / core::f64::consts::PI;
    RocketFlightForcing {
        propulsion: if powered {
            RocketFlightPropulsion::Powered {
                nozzle: RocketNozzle {
                    throat_area_m2: 1e-6,
                    exit_area_m2: exit_area,
                },
                thrust_direction_body: Vec3::new(1.0, 0.0, 0.0),
                thrust_application_arm_body_m: Vec3::new(-0.1, 0.0, 0.0),
                // Uniform circular outlet on the principal x axis. This is a
                // separately declared outlet spin moment, not tank inertia.
                co_rotating_exhaust_specific_inertia_m2: Vec3::new(
                    0.5 * exit_radius_squared,
                    0.01 + 0.25 * exit_radius_squared,
                    0.01 + 0.25 * exit_radius_squared,
                ),
                exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
            }
        } else {
            RocketFlightPropulsion::ClosedValve
        },
        additional_force_world_n: Vec3::ZERO,
        additional_torque_body_n_m: Vec3::ZERO,
    }
}

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn simulate(
    mut state: GasRigidRocketState,
    steps: usize,
    dt_s: f64,
    powered_steps: usize,
    mut output: impl Write,
) -> Result<Journey> {
    if !(1..=10000).contains(&steps)
        || powered_steps > steps
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !(steps as f64 * dt_s).is_finite()
    {
        return Err("require 1..10000 steps, 0..steps powered steps and positive finite dt".into());
    }
    let air = atmosphere()?;
    let initial_mass = state.reservoir().mass_kg();
    let (mut impulse, mut expelled) = (Sum::default(), Sum::default());
    let (mut mass_error, mut impulse_error, mut angular_error): (f64, f64, f64) = (0.0, 0.0, 0.0);
    let mut coast_mass = None;
    writeln!(
        output,
        "time_s,load_time_s,mode,nozzle_regime,altitude_m,ambient_pressure_pa,mach,reynolds,dynamic_pressure_pa,gas_mass_kg,gas_temperature_k,gas_pressure_pa,radial_velocity_m_s,attitude_w,attitude_x,attitude_y,attitude_z,angular_speed_rad_s,thrust_n,total_axial_impulse_n_s,expelled_mass_kg,total_mass_residual_kg,step_mass_residual_kg,step_axial_impulse_residual_n_s,step_angular_momentum_residual_kg_m2_s"
    )?;
    for i in 0..steps {
        let powered = i < powered_steps;
        if i == powered_steps {
            coast_mass = Some(state.reservoir().mass_kg());
        }
        let candidate = advance_atmospheric_rocket(
            &state,
            &air,
            geometry(),
            forcing(powered),
            dt_s,
            controls(),
            // Explicit manufactured drag law. A real consumer supplies its
            // validated Cd(Mach, Reynolds, body relative-flow direction, ...).
            |_| Ok(DRAG_COEFFICIENT),
            || false,
        )
        .map_err(|error| {
            format!(
                "flight interval {i} at t={} s refused: {error}",
                i as f64 * dt_s
            )
        })?;
        let (next, step, midpoint, endpoint) = candidate.into_parts();
        let mechanics = step.mechanics();
        let angular_residual = norm(mechanics.angular_momentum_residual_world_kg_m2_s);
        angular_error = angular_error.max(angular_residual);
        let (mode, regime, thrust, step_mass_error, step_impulse_error) = match &step {
            RocketFlightStep::Powered(report) => {
                impulse.add(report.mechanics.axial_thrust_impulse_n_s);
                expelled.add(report.gas.expelled_mass_kg);
                let residual = report
                    .remaining_mass_residual_kg
                    .abs()
                    .max(report.consumed_mass_residual_kg.abs());
                (
                    "powered",
                    format!("{:?}", report.midpoint_nozzle.regime),
                    report.midpoint_nozzle.thrust_n,
                    residual,
                    report.axial_impulse_residual_n_s,
                )
            }
            RocketFlightStep::Coast(_) => ("closed_valve", "none".to_owned(), 0.0, 0.0, 0.0),
        };
        mass_error = mass_error.max(step_mass_error);
        impulse_error = impulse_error.max(step_impulse_error.abs());
        state = next;
        let motion = state.vehicle().translation();
        let gas = state.reservoir().state().state();
        let radial_velocity = motion
            .velocity_m_s
            .dot(motion.position_m.scale(norm(motion.position_m).recip()));
        let q = state.vehicle().orientation().components();
        let omega = state
            .vehicle()
            .mass_properties()?
            .angular_velocity_body_checked(state.vehicle().angular_momentum_body_kg_m2_s())?;
        writeln!(
            output,
            "{:.9},{:.9},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.9},{:.9},{:.12e},{:.12e},{:.3e},{:.3e},{:.3e},{:.3e}",
            (i + 1) as f64 * dt_s,
            (i as f64 + 0.5) * dt_s,
            mode,
            regime,
            endpoint.center_altitude_m,
            midpoint.ambient_pressure_pa,
            midpoint.airflow.mach_number(),
            midpoint.airflow.reynolds_number(),
            midpoint.airflow.dynamic_pressure_pa(),
            state.reservoir().mass_kg(),
            gas.temperature,
            gas.pressure,
            radial_velocity,
            q[0],
            q[1],
            q[2],
            q[3],
            norm(omega),
            thrust,
            impulse.value,
            expelled.value,
            initial_mass - state.reservoir().mass_kg() - expelled.value,
            step_mass_error,
            step_impulse_error,
            angular_residual
        )?;
    }
    Ok(Journey {
        mass_residual_kg: initial_mass - state.reservoir().mass_kg() - expelled.value,
        state,
        axial_impulse_n_s: impulse.value,
        expelled_mass_kg: expelled.value,
        maximum_mass_residual_kg: mass_error,
        maximum_axial_impulse_residual_n_s: impulse_error,
        maximum_angular_momentum_residual_kg_m2_s: angular_error,
        coast_start_mass_kg: coast_mass,
    })
}

fn main() -> Result<()> {
    let words: Vec<_> = std::env::args().skip(1).take(7).collect();
    if words == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !(3..=6).contains(&words.len()) {
        return Err(USAGE.into());
    }
    let steps = words.get(3).map_or(Ok(200), |x| x.parse::<usize>())?;
    let dt = words.get(4).map_or(Ok(0.01), |x| x.parse::<f64>())?;
    let powered_steps = words.get(5).map_or(Ok(steps / 2), |x| x.parse::<usize>())?;
    let mut text = String::new();
    std::fs::File::open(&words[0])?
        .take(65537)
        .read_to_string(&mut text)?;
    if text.len() > 65536 {
        return Err("synthetic gas source exceeds 64 KiB".into());
    }
    let gas = ClosedGasParcel::new(
        source::load(&text, &words[0])?,
        source::point(words[1].parse()?, words[2].parse()?)?,
        TANK_VOLUME_M3,
    )?;
    let state = initial(gas)?;
    eprintln!(
        "Manufactured 2 L gas tank, 0.1 kg dry body, 1 mm2 throat, area ratio 2; initially airborne at 100 m with 30 m/s radial velocity. Explicit Cd=0.5, spherical standard atmosphere, 24-hour rotation, 5 m/s background wind. Valve closes after {powered_steps} steps. This is a model-composition example, not calibrated vehicle performance."
    );
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    let result = simulate(state, steps, dt, powered_steps, &mut output)?;
    output.flush()?;
    let altitude = norm(result.state.vehicle().translation().position_m) - USSA_EFFECTIVE_RADIUS_M;
    eprintln!(
        "Completed {:.6} s; altitude={:.9} m, axial impulse={:.12e} N s, expelled={:.12e} kg, mass residual={:.3e} kg. Maximum step residuals: mass={:.3e} kg, axial impulse={:.3e} N s, angular momentum={:.3e} kg m2/s. Coast start mass={:?} kg.",
        steps as f64 * dt,
        altitude,
        result.axial_impulse_n_s,
        result.expelled_mass_kg,
        result.mass_residual_kg,
        result.maximum_mass_residual_kg,
        result.maximum_axial_impulse_residual_n_s,
        result.maximum_angular_momentum_residual_kg_m2_s,
        result.coast_start_mass_kg
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complete_source_nozzle_attitude_atmosphere_and_valve_journey() {
        let input = include_str!("../../fs-material/examples/material_gas_spring/mixture-a.tsv");
        let gas = ClosedGasParcel::new(
            source::load(input, "manufactured atmospheric journey").unwrap(),
            source::point(1000.0, 1e6).unwrap(),
            TANK_VOLUME_M3,
        )
        .unwrap();
        let start = initial(gas).unwrap();
        let initial_orientation = start.vehicle().orientation();
        let result = simulate(start, 20, 0.01, 10, std::io::sink()).unwrap();
        assert!(result.axial_impulse_n_s > 0.0 && result.expelled_mass_kg > 0.0);
        assert!(result.mass_residual_kg.abs() < 1e-13);
        assert!(result.maximum_mass_residual_kg < 1e-13);
        assert!(result.maximum_axial_impulse_residual_n_s < 1e-14);
        assert!(result.maximum_angular_momentum_residual_kg_m2_s < 1e-15);
        assert_eq!(
            result.coast_start_mass_kg.unwrap(),
            result.state.reservoir().mass_kg()
        );
        assert_ne!(initial_orientation, result.state.vehicle().orientation());
        assert!(
            norm(result.state.vehicle().translation().position_m) > USSA_EFFECTIVE_RADIUS_M + 100.0
        );
    }
}
