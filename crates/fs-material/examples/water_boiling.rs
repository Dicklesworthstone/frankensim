//! Boil and condense water with the shared thermal timeline and IAPWS source.
//!
//! Run: `cargo run -p fs-material --example water_boiling`.
//! The water starts on its saturation plateau at atmospheric pressure. A
//! heater runs for 600 s, then prescribed cooling reverses the net heat while
//! ambient heat loss continues. The external pressure load permits expansion;
//! the cell geometry is a fixed reference heat-exchange area, not a rigid tank.
//! This example stays inside liquid/vapor coexistence throughout. Subcooling,
//! superheating, interface kinetics and changing contact area need other laws.

use fs_material::phase::vaporization::{IsobaricWaterVaporization, WaterVaporizationState};
use fs_qty::Pressure;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use fs_thermochem::water::{SOURCE_RELEASE, SOURCE_URL, WaterSaturationSolve};

const MASS_KG: f64 = 0.5;
const INITIAL_VAPOR_FRACTION: f64 = 0.2;
const HEATER_W: f64 = 1000.0;
const HEATING_END_S: f64 = 600.0;
const AMBIENT_K: f64 = 298.15;
const AREA_M2: f64 = 0.25;
const HTC_W_M2_K: f64 = 8.0;
const ENERGY_TOLERANCE_J: f64 = 1e-6;

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: AMBIENT_K,
        htc_w_m2_k: HTC_W_M2_K,
        reservoir: None,
        source_w,
    }
}

fn row(time_s: f64, state: WaterVaporizationState, net_heat_j: f64, work_j: f64, delta_u_j: f64) {
    println!(
        "{time_s:.0},{:.6},{:.9},{:.9},{:.9},{:.9},{net_heat_j:.6},{work_j:.6},{delta_u_j:.6}",
        state.temperature_k(),
        state.vapor_mass_fraction(),
        MASS_KG * state.liquid_mass_fraction(),
        MASS_KG * state.vapor_mass_fraction(),
        MASS_KG * state.specific_volume_m3_kg(),
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = IsobaricWaterVaporization::at_pressure(
        Pressure::new(101_325.0),
        WaterSaturationSolve {
            pressure_tolerance_pa: 1e-6,
            maximum_iterations: 96,
        },
        || false,
    )?;
    let source = law.saturation();
    let initial_h = source.liquid_specific_enthalpy().value()
        + INITIAL_VAPOR_FRACTION
            * (source.vapor_specific_enthalpy().value()
                - source.liquid_specific_enthalpy().value());
    let initial = law.state_at_specific_enthalpy(initial_h)?;
    let loss_w = HTC_W_M2_K * AREA_M2 * (initial.temperature_k() - AMBIENT_K);
    let events = vec![
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
                    surface_area_m2: AREA_M2,
                },
                boundary: boundary(HEATER_W),
            }),
        },
        TimedEvent {
            time_s: HEATING_END_S,
            // Net heat becomes exactly the negative of the heating-stage
            // power, accounting for continuing heat loss to the ambient.
            event: Event::Boundary(CellId(1), boundary(-HEATER_W + 2.0 * loss_w)),
        },
    ];
    let simulation = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 1000,
        },
    )?;
    let step = StepConfig {
        max_dt_s: 2.0,
        energy_tolerance_j: ENERGY_TOLERANCE_J,
    };
    let (mut simulation, _) = simulation.advance_to(0.0, step, || Ok(()))?;
    let mut net_heat = 0.0;
    let mut entropy_generation = 0.0;
    let mut maximum_energy_residual = 0.0_f64;
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!(
        "# p={:.9} Pa; latent_h={:.6} J/kg; ambient_loss={loss_w:.6} W",
        initial.pressure_pa(),
        source.latent_enthalpy().value(),
    );
    println!(
        "time_s,temperature_K,vapor_mass_fraction,liquid_kg,vapor_kg,volume_m3,net_heat_J,expansion_work_J,internal_energy_change_J"
    );
    row(0.0, initial, 0.0, 0.0, 0.0);
    for index in 1..=12 {
        let time = 100.0 * f64::from(index);
        let (next, heat, entropy) = simulation.advance_to_with_entropy(time, step, || Ok(()))?;
        net_heat += heat.external_heat_j;
        entropy_generation += entropy.conduction_generation_j_k;
        let actual_h = next.cells()[&CellId(1)].specific_h_j_kg;
        let actual = law.state_at_specific_enthalpy(actual_h)?;
        // The original material owner independently accounts for real phase
        // inventory, internal energy and pressure work for the supplied heat.
        let calorimetry =
            law.advance_heat(initial, MASS_KG, net_heat, ENERGY_TOLERANCE_J, || false)?;
        let trajectory_residual = MASS_KG * (actual_h - initial_h) - net_heat;
        let owner_disagreement = MASS_KG * (actual_h - calorimetry.state.specific_enthalpy_j_kg());
        for residual in [
            heat.residual_j,
            trajectory_residual,
            owner_disagreement,
            calorimetry.heat_residual_j,
            calorimetry.first_law_residual_j,
        ] {
            maximum_energy_residual = maximum_energy_residual.max(residual.abs());
            if residual.abs() > ENERGY_TOLERANCE_J {
                return Err(format!("water trajectory energy residual {residual} J").into());
            }
        }
        row(
            time,
            actual,
            net_heat,
            calorimetry.boundary_work_j,
            calorimetry.internal_energy_change_j,
        );
        simulation = next;
    }
    if net_heat.abs() > ENERGY_TOLERANCE_J || entropy_generation <= 0.0 {
        return Err("water cycle did not close or ambient heat flow violated entropy".into());
    }
    println!(
        "# max_energy_residual={maximum_energy_residual:.9e} J; ambient_conduction_entropy={entropy_generation:.9} J/K"
    );
    println!("# The prescribed heater/cooler does not specify its reservoir entropy.");
    Ok(())
}
