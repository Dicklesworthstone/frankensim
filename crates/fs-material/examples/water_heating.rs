//! Heat liquid water through boiling into superheated steam, then condense it.
//!
//! Run: `cargo run -p fs-material --example water_heating`.
//! Original IAPWS IF97 Regions 1, 2 and 4 supply all caloric properties:
//! <https://iapws.org/technical-guidance/release/IF97-Rev>.
//! A fixed mass expands against a constant atmospheric pressure load. The
//! shared thermal timeline applies prescribed heater/cooler power; physical
//! volume and pressure work come from the same material owner. The timeline's
//! geometry is reference heat-exchange geometry, not a rigid sealed vessel.
//! The prescribed source has no specified reservoir entropy. The printed
//! pressure-specific thermal entropy includes the separately reported IF97
//! saturation matching correction; it is not cross-pressure absolute entropy.

use fs_material::phase::water::{IsobaricWater, WaterChartControls, WaterState};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

const MASS_KG: f64 = 0.05;
const HEATING_END_S: f64 = 100.0;
const ENERGY_TOLERANCE_J: f64 = 1e-6;

fn boundary(power_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power_w,
    }
}

fn row(time_s: f64, state: WaterState, heat_j: f64, work_j: f64, internal_j: f64) {
    println!(
        "{time_s:.0},{:?},{:.6},{:.9},{:.9},{:.9},{heat_j:.6},{work_j:.6},{internal_j:.6}",
        state.phase(),
        state.temperature_k(),
        state.vapor_mass_fraction(),
        MASS_KG * state.vapor_mass_fraction(),
        MASS_KG * state.specific_volume_m3_kg(),
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = IsobaricWater::new(
        Pressure::new(101_325.0),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 1073.15,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65535,
            },
        },
        || false,
    )?;
    let initial = law.state_at_temperature(Temperature::new(300.0), || false)?;
    let hot = law.state_at_temperature(Temperature::new(700.0), || false)?;
    let initial_h = initial.specific_enthalpy_j_kg();
    let heating_j = MASS_KG * (hot.specific_enthalpy_j_kg() - initial_h);
    let power_w = heating_j / HEATING_END_S;
    let simulation = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(1),
                    segment_id: 1,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: MASS_KG,
                    specific_h_j_kg: initial_h,
                    geometry: CellGeometry {
                        start_m: [0.0, 0.0, 0.0],
                        end_m: [0.0, 0.0, 1.0],
                        volume_m3: MASS_KG * initial.specific_volume_m3_kg(),
                        surface_area_m2: 1.0,
                    },
                    boundary: boundary(power_w),
                }),
            },
            TimedEvent {
                time_s: HEATING_END_S,
                event: Event::Boundary(CellId(1), boundary(-power_w)),
            },
        ],
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 512,
        },
    )?;
    let config = StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: ENERGY_TOLERANCE_J,
    };
    let (mut simulation, _) = simulation.advance_to(0.0, config, || Ok(()))?;
    println!("# Original IAPWS IF97 (2012): https://iapws.org/technical-guidance/release/IF97-Rev");
    println!(
        "# mass={MASS_KG} kg; p={} Pa; saturation={:.9} K; heater={power_w:.9} W",
        initial.pressure_pa(),
        law.saturation().temperature().value(),
    );
    println!(
        "# source_entropy_seam={:.12e} J/(kg K); thermal_S=raw_s+vapor_fraction*seam",
        law.saturation_entropy_defect_j_kg_k(),
    );
    println!(
        "time_s,phase,temperature_K,vapor_mass_fraction,vapor_kg,physical_volume_m3,cumulative_heat_J,expansion_work_J,internal_energy_change_J"
    );
    row(0.0, initial, 0.0, 0.0, 0.0);
    let mut net_heat = 0.0;
    let mut thermal_entropy = 0.0;
    let mut material_entropy = 0.0;
    let mut expansion_work = 0.0;
    let mut internal_energy_change = 0.0;
    let mut max_residual = 0.0_f64;
    for time in [
        5.0, 20.0, 50.0, 80.0, 90.0, 100.0, 110.0, 120.0, 150.0, 180.0, 195.0, 200.0,
    ] {
        let (next, energy, entropy) =
            simulation.advance_to_with_entropy(time, config, || Ok(()))?;
        net_heat += energy.external_heat_j;
        thermal_entropy += entropy.thermal_storage_change_j_k;
        let h = next.cells()[&CellId(1)].specific_h_j_kg;
        let state = law.state_at_specific_enthalpy(h, || false)?.into_state();
        let previous = law
            .state_at_specific_enthalpy(simulation.cells()[&CellId(1)].specific_h_j_kg, || false)?
            .into_state();
        // Audit each finite accepted interval. A tiny signed roundoff in the
        // sum of a closed cycle's heat is not another physical heating request.
        let balance = law.advance_heat(
            previous,
            MASS_KG,
            energy.external_heat_j,
            ENERGY_TOLERANCE_J,
            || false,
        )?;
        expansion_work += balance.boundary_work_j;
        internal_energy_change += balance.internal_energy_change_j;
        material_entropy += balance.thermal_entropy_change_j_k;
        for residual in [
            energy.residual_j,
            MASS_KG * (h - initial_h) - net_heat,
            MASS_KG * (h - balance.state.specific_enthalpy_j_kg()),
            internal_energy_change + expansion_work - net_heat,
            balance.enthalpy_residual_j,
            balance.energy_residual_j,
            balance.represented_enthalpy_residual_j,
            balance.represented_energy_residual_j,
        ] {
            max_residual = max_residual.max(residual.abs());
            if residual.abs() > ENERGY_TOLERANCE_J {
                return Err(
                    format!("heat/pressure-work balance exceeds budget: {residual} J").into(),
                );
            }
        }
        if (thermal_entropy - material_entropy).abs() > 1e-8 {
            return Err("thermal entropy does not match the original material path".into());
        }
        row(
            time,
            state,
            net_heat,
            expansion_work,
            internal_energy_change,
        );
        simulation = next;
    }
    let final_h = simulation.cells()[&CellId(1)].specific_h_j_kg;
    let final_state = law
        .state_at_specific_enthalpy(final_h, || false)?
        .into_state();
    if net_heat.abs() > ENERGY_TOLERANCE_J || (final_state.temperature_k() - 300.0).abs() > 1e-7 {
        return Err("the heating/condensation cycle did not close".into());
    }
    println!(
        "# heating={heating_j:.9} J; max_energy_residual={max_residual:.9e} J; cycle_heat={net_heat:.9e} J; cycle_thermal_entropy={thermal_entropy:.9e} J/K"
    );
    println!("# Positive pressure work is outward. Delta H = heat; Delta U + p Delta V = heat.");
    println!("# Source power is prescribed; a global reservoir entropy balance is not specified.");
    Ok(())
}
