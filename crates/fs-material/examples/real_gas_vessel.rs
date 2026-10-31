//! Real-gas heating and reversible volume work with explicit synthetic sources.
//!
//! Run: cargo run -p fs-material --example real_gas_vessel -- 5000000 10000 10 0.8
//! Arguments are initial pressure [Pa], total signed heat [J], steps per
//! operation, and optional final/initial volume ratio (default 1: heat only).
//! The 10-litre vessel starts at 300 K. This manufactured mixture demonstrates
//! the conservative coupling; it is not an empirical gas or vessel safety model.

use fs_material::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::ClosedRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
// Reuse the existing example's source-card loader and its bounded input format.
#[path = "material_gas_spring/source.rs"]
mod source;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4 {
        return Err("usage: real_gas_vessel [initial_pressure_pa] [total_heat_j] [steps] [final_volume_ratio]".into());
    }
    let pressure: f64 = args.first().map_or(Ok(5.0e6), |value| value.parse())?;
    let heat: f64 = args.get(1).map_or(Ok(10_000.0), |value| value.parse())?;
    let steps: usize = args.get(2).map_or(Ok(10), |value| value.parse())?;
    let volume_ratio: f64 = args.get(3).map_or(Ok(1.0), |value| value.parse())?;
    if !heat.is_finite()
        || !(1..=100_000).contains(&steps)
        || !volume_ratio.is_finite()
        || volume_ratio <= 0.0
    {
        return Err(
            "finite signed heat, positive finite volume ratio and 1..=100000 steps required".into(),
        );
    }
    let gas = source::load(
        "A 0.75 0.028 3.5 0.001 1.7e-5 0.025\nB 0.25 0.004 2.5 0 2e-5 0.15\n",
        "real_gas_vessel manufactured binary mixture",
    )?;
    let model = ThermochemicalPengRobinsonModel::new(
        gas,
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
    let mut vessel = ClosedRealGasParcel::new(
        model,
        source::point(300.0, pressure)?,
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )?;
    let controls = GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 1000.0,
        energy_tolerance_j_kg: 1.0e-7,
        maximum_iterations: 128,
    };
    let entropy_controls = GasEntropySolve {
        minimum_temperature_k: controls.minimum_temperature_k,
        maximum_temperature_k: controls.maximum_temperature_k,
        entropy_tolerance_j_kg_k: 1.0e-9,
        maximum_iterations: controls.maximum_iterations,
    };
    println!(
        "step,supplied_heat_j,temperature_k,pressure_pa,density_kg_m3,mass_kg,cv_j_kg_k,sound_speed_m_s,energy_residual_j,volume_m3,work_by_gas_j,entropy_residual_j_k,operation"
    );
    let mut supplied = 0.0;
    let mut residual = 0.0;
    for step in 0..=steps {
        if step > 0 {
            let report = vessel.add_heat(heat / steps as f64, controls, || false)?;
            supplied += report.supplied_heat_j;
            residual += report.residual_j;
        }
        emit(&vessel, step, supplied, residual, 0.0, 0.0, "heat");
    }
    if volume_ratio != 1.0 {
        let initial_volume = vessel.volume_m3();
        let mut work = 0.0;
        let mut entropy_residual = 0.0;
        for step in 1..=steps {
            let volume = initial_volume * (1.0 + (volume_ratio - 1.0) * step as f64 / steps as f64);
            let report = vessel.set_volume_isentropic(volume, entropy_controls, || false)?;
            work += report.work_by_gas_j;
            entropy_residual += report.entropy_change_j_k;
            emit(
                &vessel,
                steps + step,
                supplied,
                residual,
                work,
                entropy_residual,
                "isentropic_volume",
            );
        }
    }
    Ok(())
}

fn emit(
    vessel: &ClosedRealGasParcel,
    step: usize,
    supplied: f64,
    residual: f64,
    work: f64,
    entropy_residual: f64,
    operation: &str,
) {
    let state = vessel.state();
    println!(
        "{step},{supplied:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{residual:.12e},{:.12e},{work:.12e},{entropy_residual:.12e},{operation}",
        state.temperature_k(),
        state.pressure_pa(),
        state.density_kg_m3(),
        vessel.mass_kg(),
        state.specific_heat_cv_j_kg_k(),
        state.sound_speed_m_s(),
        vessel.volume_m3(),
    );
}
