//! Fixed-volume real-gas heating with explicit synthetic source cards.
//!
//! Run: cargo run -p fs-material --example real_gas_vessel -- 5000000 10000 10
//! Arguments are initial pressure [Pa], total signed heat [J], and heat steps.
//! The 10-litre vessel starts at 300 K. This manufactured mixture demonstrates
//! the conservative coupling; it is not an empirical gas or vessel safety model.

use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::IsochoricRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
// Reuse the existing example's source-card loader and its bounded input format.
#[path = "material_gas_spring/source.rs"]
mod source;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 3 {
        return Err("usage: real_gas_vessel [initial_pressure_pa] [total_heat_j] [steps]".into());
    }
    let pressure: f64 = args.first().map_or(Ok(5.0e6), |value| value.parse())?;
    let heat: f64 = args.get(1).map_or(Ok(10_000.0), |value| value.parse())?;
    let steps: usize = args.get(2).map_or(Ok(10), |value| value.parse())?;
    if !heat.is_finite() || !(1..=100_000).contains(&steps) {
        return Err("finite signed heat and 1..=100000 steps required".into());
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
    let mut vessel = IsochoricRealGasParcel::new(
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
    println!(
        "step,supplied_heat_j,temperature_k,pressure_pa,density_kg_m3,mass_kg,cv_j_kg_k,sound_speed_m_s,energy_residual_j"
    );
    let mut supplied = 0.0;
    let mut residual = 0.0;
    for step in 0..=steps {
        if step > 0 {
            let report = vessel.add_heat(heat / steps as f64, controls, || false)?;
            supplied += report.supplied_heat_j;
            residual += report.residual_j;
        }
        let state = vessel.state();
        println!(
            "{step},{supplied:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{residual:.12e}",
            state.temperature_k(),
            state.pressure_pa(),
            state.density_kg_m3(),
            vessel.mass_kg(),
            state.specific_heat_cv_j_kg_k(),
            state.sound_speed_m_s()
        );
    }
    Ok(())
}
