//! Executable two-pass/layered-strip reduced thermal example, SI units.
//!
//! Run with `cargo run -p fs-conduction --example deposited_strip --release`.
//! The path import also permits `rustc --edition 2024 -O
//! crates/fs-conduction/examples/deposited_strip.rs -o deposited_strip` when
//! the workspace's external constellation is unavailable. It uses the exact
//! production kernel, with no copied implementation or substitute dependencies.
//! Coefficients are declared illustrative inputs, not calibrated filament data.

#[path = "../src/transient/deposition.rs"]
#[allow(dead_code)]
mod deposition;

use deposition::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, ThermalLink,
    TimedEvent,
};

struct DeclaredSensiblePolymer;
impl ThermalLaw for DeclaredSensiblePolymer {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(-200_000.0..=600_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h / 2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 2000.0
    }
}

fn at(time_s: f64, event: Event) -> TimedEvent {
    TimedEvent { time_s, event }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = DeclaredSensiblePolymer;
    let mut schedule = Vec::new();
    // Four cells per layer; four layers. First/second layers are the two passes.
    // Bodies and segments belong to the consumer, not the thermal solver.
    for layer in 0..4_u64 {
        for x in 0..4_u64 {
            let id = layer * 4 + x;
            let time_s = layer as f64 * 2.0 + x as f64 * 0.25;
            let boundary = Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 20.0,
                reservoir: if layer == 0 {
                    Some(Reservoir {
                        area_m2: 4e-6,
                        conductance_w_k: 0.004,
                        temperature_k: 330.0,
                    })
                } else {
                    None
                },
                source_w: 0.0,
            };
            schedule.push(at(
                time_s,
                Event::Birth(DepositedCell {
                    id: CellId(id),
                    segment_id: id,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: 8e-7,
                    specific_h_j_kg: 400_000.0,
                    geometry: CellGeometry {
                        start_m: [x as f64 * 0.01, 0.0, layer as f64 * 0.0002],
                        end_m: [(x + 1) as f64 * 0.01, 0.0, layer as f64 * 0.0002],
                        volume_m3: 8e-10,
                        surface_area_m2: 12.16e-6,
                    },
                    boundary,
                }),
            ));
            if layer > 0 {
                schedule.push(at(
                    time_s,
                    Event::Contact(ThermalLink {
                        id: 100 + id,
                        a: CellId(id - 4),
                        b: CellId(id),
                        area_m2: 4e-6,
                        conductance_w_k: 0.004,
                    }),
                ));
            }
            if x > 0 {
                schedule.push(at(
                    time_s,
                    Event::Contact(ThermalLink {
                        id: 200 + id,
                        a: CellId(id - 1),
                        b: CellId(id),
                        area_m2: 8e-8,
                        conductance_w_k: 1.6e-6,
                    }),
                ));
            }
        }
    }
    // Prescribed bed release is a boundary input, not predicted detachment.
    for id in 0..4 {
        schedule.push(at(
            15.0,
            Event::Boundary(
                CellId(id),
                Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 20.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            ),
        ));
    }
    let mut sim = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 16,
            max_links: 32,
            max_events: 64,
            max_substeps: 100_000,
        },
    )?;
    let config = StepConfig {
        max_dt_s: 0.002,
        energy_tolerance_j: 1e-10,
    };
    let mut incoming_mass = 0.0;
    let mut incoming_energy = 0.0;
    let mut external_heat = 0.0;
    let mut first_before = 0.0;
    let mut first_reheated = 0.0_f64;
    println!(
        "time_s,cells,mass_kg,stored_energy_j,first_cell_K,exposed_first_m2,advance_residual_j"
    );
    for target in [0.0, 1.75, 2.25, 4.0, 6.0, 8.0, 15.0, 30.0] {
        let (next, audit) = sim.advance_to(target, config, || Ok(()))?;
        incoming_mass += audit.incoming_mass_kg;
        incoming_energy += audit.incoming_energy_j;
        external_heat += audit.external_heat_j;
        let t = next.temperature_k(CellId(0))?;
        if target == 1.75 {
            first_before = t;
        }
        if target == 2.25 {
            first_reheated = t;
        }
        println!(
            "{target:.2},{},{:.9e},{:.9},{t:.6},{:.9e},{:.3e}",
            next.cells().len(),
            next.mass_kg(),
            next.energy_j(),
            next.exposed_area_m2(CellId(0))?,
            audit.residual_j
        );
        sim = next;
    }
    assert!(
        first_reheated > first_before,
        "second pass must reheat the first"
    );
    assert!((sim.mass_kg() - incoming_mass).abs() < 1e-15);
    assert!((sim.energy_j() - incoming_energy - external_heat).abs() < 1e-10);
    assert_eq!(sim.cells().len(), 16);
    eprintln!(
        "Verified 16 births; reheating {:.3} -> {:.3} K; mass={:.9e} kg; energy residual={:.3e} J",
        first_before,
        first_reheated,
        incoming_mass,
        sim.energy_j() - incoming_energy - external_heat
    );
    Ok(())
}
