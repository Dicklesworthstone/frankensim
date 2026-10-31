//! A heated real-gas chamber driving the native inertial piston.
//!
//! Run: cargo run -p fs-couple --example real_gas_piston -- 5000000 100 200
//! Optional arguments are initial pressure [Pa], total signed heat [J], and
//! steps. Heat is supplied uniformly over the first quarter of the steps;
//! the remaining steps show the continuing gas/wall motion with no added heat.
//! The 1-litre chamber starts at 300 K with a stationary 2-kg, 0.01-m2 wall.
//! Exterior pressure stays at its initial value and each timestep is 50 us.
//!
//! The binary caloric, transport and Peng-Robinson inputs are synthetic. This
//! demonstrates conservative coupling, not a calibrated gas or actuator model.
//! The original source and acoustic-clock checks admit every proposed step.

use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::IsochoricRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

// Reuse the existing seven-column synthetic source-card format and loader.
#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

const INITIAL_TEMPERATURE_K: f64 = 300.0;
const INITIAL_VOLUME_M3: f64 = 0.001;
const WALL_MASS_KG: f64 = 2.0;
const AREA_M2: f64 = 0.01;
const TIMESTEP_S: f64 = 0.00005;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 3 {
        return Err("usage: real_gas_piston [initial_pressure_pa] [total_heat_j] [steps]".into());
    }
    let pressure: f64 = args.first().map_or(Ok(5.0e6), |value| value.parse())?;
    let heat: f64 = args.get(1).map_or(Ok(100.0), |value| value.parse())?;
    let steps: usize = args.get(2).map_or(Ok(200), |value| value.parse())?;
    if !pressure.is_finite()
        || pressure <= 0.0
        || !heat.is_finite()
        || !(1..=100_000).contains(&steps)
    {
        return Err(
            "positive finite pressure, finite signed heat and 1..=100000 steps required".into(),
        );
    }

    let source = source::load(
        "A 0.75 0.028 3.5 0.001 1.7e-5 0.025\nB 0.25 0.004 2.5 0 2e-5 0.15\n",
        "real_gas_piston manufactured binary mixture",
    )?;
    let model = ThermochemicalPengRobinsonModel::new(
        source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A")?,
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B")?,
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )?;
    let gas = IsochoricRealGasParcel::new(
        model,
        source::point(INITIAL_TEMPERATURE_K, pressure)?,
        INITIAL_VOLUME_M3,
        PengRobinsonPhaseV1::Vapor,
    )?;
    let initial_density = gas.density_kg_m3();
    let mut wall = GasPiston::new(gas, WALL_MASS_KG, AREA_M2, 0.0)?;
    let energy = GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 500.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let heating_steps = (steps / 4).max(1);
    let mut supplied_heat = 0.0;
    let mut gas_work = 0.0;
    println!(
        "step,time_s,temperature_k,pressure_pa,density_kg_m3,volume_m3,wall_velocity_m_s,gas_mass_kg,sound_speed_m_s,supplied_heat_j,gas_work_j,gas_internal_energy_change_j,wall_kinetic_energy_change_j,exterior_pressure_work_j,total_energy_residual_j,step_energy_residual_j"
    );
    for step in 0..=steps {
        let step_residual = if step == 0 {
            0.0
        } else {
            let heat_j = if step < heating_steps {
                heat / heating_steps as f64
            } else if step == heating_steps {
                // Retain the requested total despite roundoff in subdivision.
                heat - supplied_heat
            } else {
                0.0
            };
            let volume = wall.gas().volume_m3();
            let report = wall.step(
                PistonInput {
                    duration_s: TIMESTEP_S,
                    heat_j,
                    exterior_pressure_pa: pressure,
                    external_force_n: 0.0,
                },
                PistonSolve {
                    volume_bracket_m3: [0.95 * volume, 1.05 * volume],
                    gas: energy,
                    momentum_tolerance_n_s: 1e-10,
                    energy_tolerance_j: 1e-7,
                    maximum_iterations: 128,
                },
                || false,
            )?;
            supplied_heat += report.gas.supplied_heat_j;
            gas_work += report.gas.work_by_gas_j;
            report.total_energy_residual_j
        };

        let gas = wall.gas();
        let state = gas.state();
        // Audit from the original state, including density-dependent PR
        // departures. Subtracting absolute formation energies can lose work.
        let internal_energy_change = gas.mass_kg()
            * gas
                .model()
                .specific_internal_energy_change_between_densities(
                    initial_density,
                    INITIAL_TEMPERATURE_K,
                    gas.density_kg_m3(),
                    state.temperature_k(),
                )?;
        let kinetic_energy_change = 0.5 * wall.mass_kg() * wall.velocity_m_s().powi(2);
        let exterior_work = pressure * (gas.volume_m3() - INITIAL_VOLUME_M3);
        let total_residual =
            (internal_energy_change + kinetic_energy_change) + exterior_work - supplied_heat;
        println!(
            "{step},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{supplied_heat:.12e},{gas_work:.12e},{internal_energy_change:.12e},{kinetic_energy_change:.12e},{exterior_work:.12e},{total_residual:.12e},{step_residual:.12e}",
            step as f64 * TIMESTEP_S,
            state.temperature_k(),
            state.pressure_pa(),
            gas.density_kg_m3(),
            gas.volume_m3(),
            wall.velocity_m_s(),
            gas.mass_kg(),
            state.sound_speed_m_s(),
        );
    }
    Ok(())
}
