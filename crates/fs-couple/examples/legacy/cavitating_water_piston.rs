//! A finite water inventory creates and collapses vapor while moving a piston.
//!
//! Run: cargo run -p fs-couple --features cavitation --example cavitating_water_piston -- 1
//! The optional positive integer refines the time intervals without changing
//! the pressure schedule. This is a lumped equilibrium chamber demonstration;
//! no nucleation kinetics, water-hammer wave or wall impact is claimed.

use fs_couple::cavitation::water_cell::piston::{
    WaterPiston, WaterPistonControls, WaterPistonInput,
};
use fs_couple::cavitation::water_cell::volume::WaterVolumeControls;
use fs_couple::cavitation::water_cell::{EquilibriumWaterCell, WaterCellControls};
use fs_material::phase::water::isochoric::IsochoricWaterControls;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let refinement: usize = std::env::args().nth(1).map_or(Ok(1), |v| v.parse())?;
    if !(1..=64).contains(&refinement) {
        return Err("refinement must lie in 1..=64".into());
    }
    let mass_kg = 1e-3;
    let liquid = If97State::at_pressure_temperature(
        Pressure::new(120_000.0),
        Temperature::new(373.15),
        If97Region::Region1,
    )?;
    let cell = EquilibriumWaterCell::new(
        mass_kg * liquid.specific_volume().value(),
        mass_kg,
        mass_kg * liquid.specific_internal_energy().value(),
        0.0,
        WaterCellControls {
            equilibrium: IsochoricWaterControls {
                minimum_temperature_k: 371.15,
                maximum_temperature_k: 375.15,
                minimum_pressure_pa: 100.0,
                maximum_pressure_pa: 100e6,
                volume_tolerance_m3_kg: 1e-13,
                internal_energy_tolerance_j_kg: 1e-7,
                maximum_pressure_iterations: 96,
                maximum_temperature_iterations: 96,
            },
            maximum_relative_mass_exchange: 0.1,
            maximum_inflow_parcels: 4,
            mass_tolerance_kg: 1e-15,
            energy_tolerance_j: 1e-9,
            volume_tolerance_m3: 1e-16,
        },
        || false,
    )?;
    let mut piston = WaterPiston::new(cell, 0.1, 1e-4, 0.0)?;
    let initial_energy = piston.water().state().internal_energy_j() + piston.kinetic_energy_j();
    let dt = 2e-6 / refinement as f64;
    let mut exterior_work = 0.0;
    let mut maximum_vapor_mass: f64 = 0.0;
    let mut maximum_residual: f64 = 0.0;
    let mut first_vapor_time = None;
    let mut collapse_time = None;
    println!(
        "time_s\tvolume_m3\tpressure_pa\ttemperature_k\twall_velocity_m_s\tvapor_mass_kg\tvapor_volume_m3\tkinetic_energy_j\texterior_work_j\ttotal_energy_residual_j"
    );
    emit(&piston, exterior_work, 0.0);
    for i in 1..=300 * refinement {
        let volume = piston.water().state().volume_m3();
        let exterior_pressure_pa = if i <= 100 * refinement {
            80_000.0
        } else {
            140_000.0
        };
        let step = piston.step(
            i as f64 * dt,
            WaterPistonInput {
                heat_j: 0.0,
                exterior_pressure_pa,
                external_force_n: 0.0,
            },
            WaterPistonControls {
                water: WaterVolumeControls {
                    maximum_step_s: 2.01e-6,
                    maximum_relative_volume_change: 0.1,
                    maximum_temperature_increment_k: 2.0,
                    energy_tolerance_j: 1e-9,
                    maximum_iterations: 96,
                },
                volume_bracket_m3: [volume * 0.999, volume * 1.001],
                momentum_tolerance_n_s: 1e-11,
                energy_tolerance_j: 1e-8,
                maximum_iterations: 96,
            },
            || false,
        )?;
        exterior_work += step.exterior_pressure_work_j;
        let state = piston.water().state();
        let vapor_mass = state.phase_masses_kg()[1];
        maximum_vapor_mass = maximum_vapor_mass.max(vapor_mass);
        if vapor_mass > 0.0 && first_vapor_time.is_none() {
            first_vapor_time = Some(state.time_s());
        }
        if vapor_mass == 0.0 && first_vapor_time.is_some() && collapse_time.is_none() {
            collapse_time = Some(state.time_s());
        }
        // Independent cumulative first law, not a sum of reported residuals.
        let residual =
            state.internal_energy_j() + piston.kinetic_energy_j() - initial_energy + exterior_work;
        maximum_residual = maximum_residual.max(residual.abs());
        if i % (5 * refinement) == 0 {
            emit(&piston, exterior_work, residual);
        }
    }
    eprintln!("water_mass_kg={mass_kg:e}; wall_mass_kg=0.1; area_m2=1e-4; dt_s={dt:e}");
    eprintln!(
        "first_vapor_time_s={first_vapor_time:?}; collapse_time_s={collapse_time:?}; maximum_vapor_mass_kg={maximum_vapor_mass:e}; maximum_cumulative_energy_residual_j={maximum_residual:e}"
    );
    Ok(())
}

fn emit(piston: &WaterPiston, exterior_work: f64, residual: f64) {
    let state = piston.water().state();
    let source = state.source_state();
    println!(
        "{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}",
        state.time_s(),
        state.volume_m3(),
        source.pressure_pa(),
        source.temperature_k(),
        piston.velocity_m_s(),
        state.phase_masses_kg()[1],
        state.phase_volumes_m3()[1],
        piston.kinetic_energy_j(),
        exterior_work,
        residual
    );
}
