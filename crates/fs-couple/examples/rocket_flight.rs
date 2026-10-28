//! Finite source gas -> nozzle -> propellant depletion -> vehicle translation.
//!
//! Run from the repository root:
//! `cargo run -p fs-couple --example rocket_flight --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv 500 300000`
//!
//! This uses the existing manufactured gas-law file, not measured propellant
//! data. The declared demonstration is a 2 L tank, 1 mm² throat, area ratio 1.5,
//! vacuum ambient, +x thrust, no gravity or other external force, and a 0.1 kg
//! dry vehicle. The source/thermodynamic domain and accepted gas impulse govern
//! every motion step. No independent thrust or consumption curve is supplied.

#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use fs_couple::rocket_blowdown::{GasRocketEnvironment, GasRocketState};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    entropy::GasEntropySolve,
    parcel::{
        ClosedGasParcel,
        transfer::nozzle::{
            blowdown::{GasBlowdownSolve, RocketBlowdownSolve},
            rocket::{RocketNozzle, RocketNozzleSolve},
        },
    },
};
use fs_mbd::Vec3;
use std::{
    error::Error,
    io::{Read, Write},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const USAGE: &str = "usage: rocket_flight SOURCE.tsv TEMPERATURE_K PRESSURE_PA [steps=100] [dt_s=0.02] [dry_mass_kg=0.1]";

fn controls() -> RocketBlowdownSolve {
    RocketBlowdownSolve {
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
    }
}

fn simulate(
    mut state: GasRocketState,
    steps: usize,
    dt_s: f64,
    mut output: impl Write,
) -> Result<GasRocketState> {
    if !(1..=10000).contains(&steps)
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !(steps as f64 * dt_s).is_finite()
    {
        return Err("require 1..10000 steps and a positive finite physical timestep".into());
    }
    let environment = GasRocketEnvironment {
        inertial_frame: 0,
        thrust_direction_world: Vec3::new(1.0, 0.0, 0.0),
        ambient_pressure_pa: 0.0,
        gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::ZERO,
    };
    let nozzle = RocketNozzle {
        throat_area_m2: 1e-6,
        exit_area_m2: 1.5e-6,
    };
    let mut flow = state
        .reservoir()
        .rocket_nozzle_flow(0.0, nozzle, controls().nozzle, || false)?;
    let (mut impulse, mut compensation) = (0.0, 0.0);
    writeln!(
        output,
        "time_s,flow_time_s,gas_mass_kg,temperature_k,pressure_pa,position_x_m,velocity_x_m_s,mass_flow_kg_s,thrust_n,total_impulse_n_s,remaining_mass_residual_kg,consumed_mass_residual_kg,impulse_residual_n_s"
    )?;
    for i in 0..=steps {
        let (mut remaining_error, mut consumed_error, mut impulse_error) = (0.0, 0.0, 0.0);
        if i > 0 {
            let candidate = state
                .advance(environment, nozzle, dt_s, controls(), || false)
                .map_err(|error| {
                    format!(
                        "joint gas/vehicle step {i} refused at t={} s: {error}",
                        (i - 1) as f64 * dt_s
                    )
                })?;
            let (next, step) = candidate.into_parts();
            let increment = step.mechanics.thrust_impulse_n_s.x - compensation;
            let total = impulse + increment;
            compensation = (total - impulse) - increment;
            impulse = total;
            remaining_error = step.remaining_mass_residual_kg;
            consumed_error = step.consumed_mass_residual_kg;
            impulse_error = step.thrust_impulse_residual_n_s.x;
            flow = step.midpoint_nozzle;
            state = next;
        }
        let gas = state.reservoir().state().state();
        let vehicle = state.vehicle();
        let flow_time = if i == 0 { 0.0 } else { (i as f64 - 0.5) * dt_s };
        writeln!(
            output,
            "{:.9},{:.9},{:.12e},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.3e},{:.3e},{:.3e}",
            i as f64 * dt_s,
            flow_time,
            state.reservoir().mass_kg(),
            gas.temperature,
            gas.pressure,
            vehicle.position_m.x,
            vehicle.velocity_m_s.x,
            flow.mass_flow_kg_s,
            flow.thrust_n,
            impulse,
            remaining_error,
            consumed_error,
            impulse_error
        )?;
    }
    Ok(state)
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
    let steps = words.get(3).map_or(Ok(100), |x| x.parse::<usize>())?;
    let dt_s = words.get(4).map_or(Ok(0.02), |x| x.parse::<f64>())?;
    let dry_mass = words.get(5).map_or(Ok(0.1), |x| x.parse::<f64>())?;
    let mut text = String::new();
    std::fs::File::open(&words[0])?
        .take(65537)
        .read_to_string(&mut text)?;
    if text.len() > 65536 {
        return Err("synthetic gas source exceeds 64 KiB".into());
    }
    let reservoir = ClosedGasParcel::new(
        source::load(&text, &words[0])?,
        source::point(words[1].parse()?, words[2].parse()?)?,
        0.002,
    )?;
    let state = GasRocketState::new(reservoir, dry_mass, Vec3::ZERO, Vec3::ZERO, 0)?;
    eprintln!(
        "Manufactured gas-law demonstration: 2 L tank, 1 mm2 throat, area ratio 1.5, vacuum, +x thrust, zero gravity/drag; dry mass {dry_mass} kg."
    );
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    let final_state = simulate(state, steps, dt_s, &mut output)?;
    output.flush()?;
    eprintln!(
        "Completed {:.6} s: x={:.9} m, vx={:.9} m/s, remaining propellant={:.12e} kg.",
        steps as f64 * dt_s,
        final_state.vehicle().position_m.x,
        final_state.vehicle().velocity_m_s.x,
        final_state.reservoir().mass_kg()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_source_nozzle_tank_and_vehicle_complete_the_same_journey() {
        let data = include_str!("../../fs-material/examples/material_gas_spring/mixture-a.tsv");
        let reservoir = ClosedGasParcel::new(
            source::load(data, "original manufactured mixture-a").unwrap(),
            source::point(500.0, 300000.0).unwrap(),
            0.002,
        )
        .unwrap();
        let initial = reservoir.mass_kg();
        let state = GasRocketState::new(reservoir, 0.1, Vec3::ZERO, Vec3::ZERO, 0).unwrap();
        let mut csv = Vec::new();
        let final_state = simulate(state, 100, 0.02, &mut csv).unwrap();
        assert_eq!(String::from_utf8(csv).unwrap().lines().count(), 102);
        assert!(final_state.vehicle().velocity_m_s.x > 0.0);
        assert!(final_state.vehicle().position_m.x > 0.0);
        assert!(final_state.reservoir().mass_kg() < initial);
        assert!(
            (final_state.reservoir().mass_kg() - final_state.vehicle().propellant_mass_kg).abs()
                < 1e-13
        );
        assert_eq!(final_state.vehicle().dry_mass_kg, 0.1);
        assert_eq!(final_state.vehicle().velocity_m_s.y, 0.0);
        assert_eq!(final_state.vehicle().velocity_m_s.z, 0.0);
    }
}
