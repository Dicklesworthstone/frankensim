//! Prescribed nuclear heat and nonlocal radiation in finite IF97 water receivers.
//!
//! Run: cargo run -p fs-nuclear --features thermal --example prescribed_fission_water
//! Nuclear coefficients and the prescribed event rate are SYNTHETIC examples,
//! without an isotope calibration. Water supplies the original caloric law;
//! placing the heat source in a water cell does not identify water as fuel.
//! The isobaric cells use reference contact geometry, not a rigid sealed vessel.
//! Radiation-response fractions are also SYNTHETIC, representing the output of
//! an external transport calculation rather than water attenuation coefficients.

use fs_material::phase::water::{IsobaricWater, WaterChartControls, WaterState};
use fs_nuclear::radiation::{
    RadiationDeposit, RadiationResponse, RadiationResponseLimits, RadiationSourceResponse,
};
use fs_nuclear::thermal::{
    FissionCellSource, FissionHeating, FissionRateChange, FissionRateSchedule, RadiationReceiver,
};
use fs_nuclear::{DecayHeatGroup, EnergyChannels, FissionBudget, FissionModel, FissionState};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, TimedEvent,
};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use std::error::Error;
use std::sync::Arc;

fn water_cell(id: u64, state: WaterState) -> DepositedCell {
    let mass_kg = 0.01;
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg,
        specific_h_j_kg: state.specific_enthalpy_j_kg(),
        geometry: CellGeometry {
            start_m: [0.02 * id as f64, 0.0, 0.0],
            end_m: [0.02 * id as f64, 0.0, 0.02],
            volume_m3: mass_kg * state.specific_volume_m3_kg(),
            surface_area_m2: 0.003,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let water = IsobaricWater::new(
        Pressure::new(101_325.0),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 500.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65_535,
            },
        },
        || false,
    )?;
    let initial = water.state_at_temperature(Temperature::new(300.0), || false)?;
    let controls = StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-6,
    };
    let scheduled = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(water_cell(1, initial)),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(water_cell(2, initial)),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.0001,
                    // Independently prescribed effective contact, not inferred
                    // from the water caloric law or a fictitious conductivity.
                    conductance_w_k: 1.0,
                }),
            },
        ],
        &[(MaterialId(1), &water)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 4_096,
        },
    )?;
    let (mut thermal, _) = scheduled.advance_to(0.0, controls, || Ok(()))?;
    let initial_thermal_j = thermal.energy_j();

    let model = Arc::new(FissionModel::new(
        "synthetic-water-heating-example-v1".into(),
        EnergyChannels {
            deposited_j: 4e-11,
            radiation_j: 1e-11,
            escaped_j: 2e-11,
        },
        vec![DecayHeatGroup {
            decay_constant_s_inv: 0.05,
            energy_per_fission: EnergyChannels {
                deposited_j: 2e-11,
                radiation_j: 1e-11,
                escaped_j: 0.0,
            },
        }],
    )?);
    let source = FissionCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        // Nuclear inventory is independent of the receiver's water mass.
        // No prior irradiation is inferred: its initial delayed energy is zero.
        state: FissionState::new(model, 1e15, vec![0.0])?,
        rates: FissionRateSchedule::new(
            1e13,
            vec![FissionRateChange {
                time_s: 20.0,
                fissions_per_s: 0.0,
            }],
        )?,
    };
    let mut heating = FissionHeating::new(
        &thermal,
        vec![source],
        FissionBudget {
            max_decay_groups: 1,
            energy_tolerance_j: 1e-8,
            fuel_tolerance_nuclei: 1.0,
        },
    )?
    .with_radiation_response(
        &thermal,
        Arc::new(RadiationResponse::new(
            "synthetic-fixed-spectrum-water-receivers-v1".into(),
            vec![RadiationSourceResponse {
                source_id: 1,
                deposits: vec![
                    RadiationDeposit {
                        receiver_id: 1,
                        fraction: 0.25,
                    },
                    RadiationDeposit {
                        receiver_id: 2,
                        fraction: 0.5,
                    },
                ],
                escape_fraction: 0.25,
            }],
            RadiationResponseLimits {
                max_sources: 1,
                max_receivers: 2,
                max_transfers: 2,
            },
        )?),
        vec![
            RadiationReceiver {
                cell_id: CellId(1),
                material_id: MaterialId(1),
            },
            RadiationReceiver {
                cell_id: CellId(2),
                material_id: MaterialId(1),
            },
        ],
    )?;
    let mut shutdown_fissions = None;
    let mut heat_after_shutdown_j = 0.0;
    let mut received_radiation_j = 0.0;
    let mut escaping_radiation_j = 0.0;

    println!(
        "time_s,receiver_1_k,receiver_2_k,receiver_1_vapor_fraction,fissions,local_deposition_j,emitted_radiation_j,received_radiation_j,escaping_radiation_j,nonrecoverable_export_j,delayed_storage_j"
    );
    for target_s in [0.0, 5.0, 20.0, 40.0, 80.0] {
        // The source and the original thermal timeline advance atomically.
        let (next_thermal, next_heating, heat) =
            heating.advance_to(&thermal, target_s, controls, || Ok(()))?;
        thermal = next_thermal;
        heating = next_heating;
        received_radiation_j += heating.audit().radiation_deposited_j;
        escaping_radiation_j += heating.audit().radiation_escaped_j;
        let source = heating
            .source(CellId(1))
            .expect("the bound source remains present");
        let fissions = source.state.cumulative_fissions();
        let released = heating.cumulative_released()?;
        if target_s == 20.0 {
            shutdown_fissions = Some(fissions);
        } else if target_s > 20.0 {
            assert_eq!(Some(fissions), shutdown_fissions);
            heat_after_shutdown_j += heat.external_heat_j;
        }
        let receiver = water
            .state_at_specific_enthalpy(thermal.cells()[&CellId(1)].specific_h_j_kg, || false)?
            .into_state();
        println!(
            "{:.1},{:.6},{:.6},{:.6},{:.6e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
            target_s,
            receiver.temperature_k(),
            thermal.temperature_k(CellId(2))?,
            receiver.vapor_mass_fraction(),
            fissions,
            released.deposited_j,
            released.radiation_j,
            received_radiation_j,
            escaping_radiation_j,
            released.escaped_j,
            source.state.stored_delayed_energy_j(),
        );
    }
    assert!(heat_after_shutdown_j > 0.0);
    assert!(thermal.temperature_k(CellId(2))? > 300.0);
    // Internal contact cancels. Radiation deposits once; escaping radiation and
    // the distinct nonrecoverable source channel never enter thermal storage.
    let released = heating.cumulative_released()?;
    let deposited_j = released.deposited_j + received_radiation_j;
    assert!((received_radiation_j + escaping_radiation_j - released.radiation_j).abs() <= 1e-6);
    assert!((thermal.energy_j() - initial_thermal_j - deposited_j).abs() <= 1e-6);
    Ok(())
}
