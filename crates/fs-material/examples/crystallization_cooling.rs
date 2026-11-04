//! Finite-rate crystallization against a finite mold on the shared thermal owner.
//!
//! Run: cargo run -p fs-material --example crystallization_cooling
//!
//! All numerical material inputs below are MANUFACTURED demonstration values,
//! not a calibrated polymer, alloy or water model. Replace BOTH original
//! calorics and the kinetic data together with applicable source-backed data.
//! Output distinguishes mass crystallinity, volume crystallinity and free
//! volume; heat-transfer geometry remains the declared reference geometry.

use fs_blake3::ContentHash;
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_thermal::crystallization::{CrystallizationRateKnot, TabulatedCrystallizationRate};
use fs_thermal::enthalpy::crystallization::{
    CrystallizationControls, CrystallizationMaterial, CrystallizationSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, TimedEvent,
};

fn caloric(
    id: u8,
    phase: SolidLiquidPhase,
    cp: f64,
    density: f64,
    reference_h: f64,
) -> Result<EquilibriumEnthalpyPhaseCurve, Box<dyn std::error::Error>> {
    Ok(EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([id; 32]),
        phase,
        reference_h,
        &[280.0, 420.0].map(|t| HeatCapacityKnot {
            temperature_k: t,
            specific_heat_capacity_j_kg_k: cp,
            bulk_density_kg_m3: density,
        }),
        1e-6,
    )?)
}

fn cell(id: u64, mass: f64, h: f64, density: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(id),
            mass_kg: mass,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [id as f64 * 0.01, 0.0, 0.0],
                end_m: [(id + 1) as f64 * 0.01, 0.0, 0.0],
                volume_m3: mass / density,
                surface_area_m2: 0.01,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                source_w: 0.0,
                reservoir: None,
            },
        }),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let parent = caloric(17, SolidLiquidPhase::Liquid, 2000.0, 950.0, 120_000.0)?;
    let crystal = caloric(17, SolidLiquidPhase::Solid, 1500.0, 1100.0, 0.0)?;
    let mold = caloric(18, SolidLiquidPhase::Solid, 900.0, 2700.0, 0.0)?;
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 420.0])?;
    let rate = TabulatedCrystallizationRate::new(
        [
            (280.0, 0.0),
            (320.0, 0.15),
            (360.0, 0.4),
            (400.0, 0.1),
            (420.0, 0.0),
        ]
        .map(|(temperature_k, rate_s_inv)| CrystallizationRateKnot {
            temperature_k,
            rate_s_inv,
        })
        .to_vec(),
    )?;
    let law = pair.law(&rate, ContentHash([19; 32]), 2.5, 0.65, 420.0)?;
    let thermal = DepositionSimulation::new(
        vec![
            cell(
                1,
                0.01,
                pair.parent().enthalpy_at_temperature(370.0)?,
                950.0,
            ),
            cell(2, 0.05, 900.0 * (300.0 - 280.0), 2700.0),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.001,
                    conductance_w_k: 0.2,
                }),
            },
        ],
        &[(MaterialId(1), pair.parent()), (MaterialId(2), &mold)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 20_000,
        },
    )?
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )?
    .0;
    let mut simulation = CrystallizationSimulation::new(
        thermal,
        &[CrystallizationMaterial {
            material_id: MaterialId(1),
            amorphous: pair.parent(),
            law: &law,
            maximum_temperature_slope: pair.maximum_temperature_slope(),
        }],
        || Ok(()),
    )?;
    let initial_energy = simulation.energy_j();
    println!(
        "time_s,material_temperature_K,mold_temperature_K,crystal_mass_fraction,crystal_volume_fraction,free_material_volume_m3,enthalpy_change_J"
    );
    for sample in 0..=20 {
        if sample > 0 {
            simulation = simulation
                .advance_to(sample as f64, CrystallizationControls::default(), || Ok(()))?
                .0;
        }
        let state = simulation.states()[&CellId(1)];
        let volume = pair.volume_state(state.temperature_k(), state.crystalline_mass_fraction())?;
        println!(
            "{:.3},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.6e}",
            simulation.time_s(),
            state.temperature_k(),
            simulation.temperature_k(CellId(2))?,
            state.crystalline_mass_fraction(),
            volume.crystalline_volume_fraction,
            simulation.cells()[&CellId(1)].mass_kg * volume.specific_volume_m3_kg,
            simulation.energy_j() - initial_energy
        );
    }
    Ok(())
}
