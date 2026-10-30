//! Freeze and thaw water with the existing enthalpy timeline and original IAPWS sources.
//!
//! Run: `cargo run -p fs-material --example water_freeze_thaw`.
//! An assumed internally equilibrated 50 g sample exchanges heat with a cold
//! reservoir for 300 s, then a warm reservoir. Conductance and reference contact
//! area are explicit illustrative inputs. The external pressure load permits
//! expansion; thermal geometry remains reference heat-transfer geometry.
//! Temperature and phase inventory follow enthalpy. This equilibrium example
//! does not resolve an ice front, nucleation, transport inside the sample, or
//! changes in contact area during freezing.

use fs_material::phase::water::freezing::thermal::{
    FreezingThermalControls, FreezingWaterThermalLaw,
};
use fs_material::phase::water::freezing::{
    FreezingWater, FreezingWaterControls, FreezingWaterState,
};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, Reservoir, StepConfig, TimedEvent,
};
use fs_thermochem::water::ice::bounds::IceIhHeatCapacityBoundSolve;
use fs_thermochem::water::ice::{
    IceMeltingSolve, MELTING_SOURCE_RELEASE, MELTING_SOURCE_URL, SOURCE_RELEASE, SOURCE_URL,
};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

const MASS_KG: f64 = 0.05;
const PRESSURE_PA: f64 = 101_325.0;
const COLD_K: f64 = 253.15;
const WARM_K: f64 = 293.15;
const CONTACT_CONDUCTANCE_W_K: f64 = 5.0;
const CONTACT_AREA_M2: f64 = 0.1;
const COOLING_END_S: f64 = 300.0;
const FINAL_TIME_S: f64 = 650.0;
const ENERGY_TOLERANCE_J: f64 = 2e-6;

fn boundary(reservoir_k: f64) -> Boundary {
    Boundary {
        ambient_k: WARM_K,
        htc_w_m2_k: 0.0,
        reservoir: Some(Reservoir {
            area_m2: CONTACT_AREA_M2,
            conductance_w_k: CONTACT_CONDUCTANCE_W_K,
            temperature_k: reservoir_k,
        }),
        source_w: 0.0,
    }
}

fn phase(state: FreezingWaterState) -> &'static str {
    match state.liquid_mass_fraction() {
        0.0 => "ice",
        1.0 => "liquid",
        _ => "ice+liquid",
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chart = FreezingWater::new(
        Pressure::new(PRESSURE_PA),
        FreezingWaterControls {
            minimum_temperature_k: 243.15,
            maximum_temperature_k: 303.15,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            maximum_gibbs_mismatch_j_kg: 30.0,
            melting: IceMeltingSolve {
                pressure_tolerance_pa: 1e-5,
                temperature_tolerance_k: 1e-10,
                maximum_iterations: 96,
            },
        },
        || false,
    )?;
    let law = FreezingWaterThermalLaw::new(
        chart,
        FreezingThermalControls {
            ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
            liquid_heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )?;
    let chart = law.chart();
    let initial = chart.state_at_temperature(Temperature::new(WARM_K), || false)?;
    let initial_h = initial.specific_enthalpy_j_kg();
    let initial_volume = MASS_KG * initial.specific_volume_m3_kg();
    let timeline = DepositionSimulation::new(
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
                        end_m: [0.1, 0.0, 0.0],
                        volume_m3: initial_volume,
                        surface_area_m2: CONTACT_AREA_M2,
                    },
                    boundary: boundary(COLD_K),
                }),
            },
            TimedEvent {
                time_s: COOLING_END_S,
                event: Event::Boundary(CellId(1), boundary(WARM_K)),
            },
        ],
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 4096,
        },
    )?;
    let config = StepConfig {
        max_dt_s: 0.5,
        energy_tolerance_j: ENERGY_TOLERANCE_J,
    };
    let (mut timeline, _) = timeline.advance_to(0.0, config, || Ok(()))?;
    let mut net_heat = 0.0;
    let mut entropy = [0.0; 4];
    let mut maximum_residual = 0.0_f64;
    let mut saw_freezing = false;
    let mut saw_cold_ice = false;
    let mut saw_thawing = false;
    let mut final_state = initial;
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!("# {MELTING_SOURCE_RELEASE}: {MELTING_SOURCE_URL}");
    println!(
        "# p={PRESSURE_PA} Pa; mass={MASS_KG} kg; G={CONTACT_CONDUCTANCE_W_K} W/K; cold={COLD_K} K; warm={WARM_K} K"
    );
    println!(
        "# melting_temperature={:.12} K; fusion_enthalpy={:.9} J/kg; source_gibbs_mismatch={:.9} J/kg; thermal_entropy_seam={:.9e} J/(kg*K)",
        chart.melting().temperature().value(),
        chart.melting().latent_enthalpy_j_kg(),
        chart.melting().gibbs_mismatch_j_kg(),
        chart.melting().entropy_defect_j_kg_k(),
    );
    println!(
        "time_s,phase,temperature_K,liquid_mass_fraction,ice_kg,liquid_kg,liquid_volume_fraction,free_volume_m3,net_heat_J,delta_u_J,pressure_work_J,thermal_entropy_J_K,reservoir_entropy_input_J_K,conductive_entropy_J_K,temporal_entropy_defect_J_K"
    );
    println!(
        "0,liquid,{WARM_K:.9},1.000000000,0.000000000,{MASS_KG:.9},1.000000000,{initial_volume:.12},0,0,0,0,0,0,0"
    );
    for index in 1..=13 {
        let time = 50.0 * f64::from(index);
        let (next, energy, interval_entropy) =
            timeline.advance_to_with_entropy(time, config, || Ok(()))?;
        net_heat += energy.external_heat_j;
        for (total, value) in entropy.iter_mut().zip([
            interval_entropy.thermal_storage_change_j_k,
            interval_entropy.reservoir_entropy_input_j_k,
            interval_entropy.conduction_generation_j_k,
            interval_entropy.temporal_defect_j_k,
        ]) {
            *total += value;
        }
        if interval_entropy.has_unclassified_sources {
            return Err("reservoir-only example acquired an unclassified source".into());
        }
        let actual_h = next.cells()[&CellId(1)].specific_h_j_kg;
        let state = chart
            .state_at_specific_enthalpy(actual_h, || false)?
            .into_state();
        let calorimetry =
            chart.advance_heat(initial, MASS_KG, net_heat, ENERGY_TOLERANCE_J, || false)?;
        for residual in [
            energy.residual_j,
            MASS_KG * (actual_h - initial_h) - net_heat,
            MASS_KG * (actual_h - calorimetry.state.specific_enthalpy_j_kg()),
            calorimetry.enthalpy_residual_j,
            calorimetry.energy_residual_j,
            calorimetry.represented_enthalpy_residual_j,
            calorimetry.represented_energy_residual_j,
        ] {
            maximum_residual = maximum_residual.max(residual.abs());
            if residual.abs() > ENERGY_TOLERANCE_J {
                return Err(format!("freezing trajectory energy residual {residual} J").into());
            }
        }
        let thermal_entropy =
            MASS_KG * chart.specific_entropy_change_j_kg_k(initial_h, actual_h, || false)?;
        if (thermal_entropy - entropy[0]).abs() > 2e-8
            || next.cells()[&CellId(1)].geometry.volume_m3 != initial_volume
        {
            return Err("freezing entropy or reference geometry changed inconsistently".into());
        }
        let x = state.liquid_mass_fraction();
        saw_freezing |= time < COOLING_END_S && x > 0.0 && x < 1.0;
        saw_cold_ice |= time <= COOLING_END_S && x == 0.0 && state.temperature_k() < 260.0;
        saw_thawing |= time > COOLING_END_S && x > 0.0 && x < 1.0;
        println!(
            "{time:.0},{},{:.9},{x:.9},{:.9},{:.9},{:.9},{:.12},{net_heat:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
            phase(state),
            state.temperature_k(),
            MASS_KG * state.ice_mass_fraction(),
            MASS_KG * x,
            state.liquid_volume_fraction(),
            MASS_KG * state.specific_volume_m3_kg(),
            calorimetry.internal_energy_change_j,
            calorimetry.boundary_work_j,
            entropy[0],
            entropy[1],
            entropy[2],
            entropy[3],
        );
        final_state = state;
        timeline = next;
    }
    if !saw_freezing
        || !saw_cold_ice
        || !saw_thawing
        || final_state.liquid_mass_fraction() != 1.0
        || final_state.temperature_k() <= 290.0
        || entropy[0] - entropy[1] <= 0.0
        || timeline.time_s() != FINAL_TIME_S
    {
        return Err("reservoir history did not complete physical freezing and thawing".into());
    }
    println!(
        "# maximum_energy_residual={maximum_residual:.9e} J; sample_plus_reservoir_entropy={:.9} J/K; temporal_entropy_defect={:.9e} J/K",
        entropy[0] - entropy[1],
        entropy[3],
    );
    println!(
        "# Entropy is the pressure-local integral dh/T; the raw source seam remains explicit."
    );
    println!(
        "# Final liquid approaches its warm reservoir in finite time; this is not an exact closed thermodynamic cycle."
    );
    Ok(())
}
