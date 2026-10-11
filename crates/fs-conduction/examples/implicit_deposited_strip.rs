//! Resolved two-material strip over a FINITE, cooling bed, using implicit thermal
//! transport. This is synthetic verification/benchmark input, NOT filament data,
//! a slicer, or a validated print predictor. All temperatures are kelvin.
//!
//! Run: cargo run --release -p fs-conduction --example implicit_deposited_strip -- 12 3 0.25 5
//! Arguments: segments/layer, layers, maximum step [s], final cooling dwell [s].
//! The JSON summary includes actual native elapsed time and operator work; CSV
//! rows contain only accepted final state. Compare dt separately from cell count.

use std::error::Error;
use std::io;
use std::time::Instant;

use fs_blake3::ContentHash;
use fs_conduction::transient::deposition as resolved;
use fs_conduction::transient::enthalpy::deposition::implicit::{
    DepositionBudget, DepositionMaterial, DepositionSimulation, DepositionStepConfig,
    EventId, MaterialId,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

fn curve(identity: u8) -> Result<EquilibriumEnthalpyPhaseCurve, Box<dyn Error>> {
    Ok(EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([identity; 32]), SolidLiquidPhase::Solid,
        [200.0, 800.0].into_iter().map(|temperature_k| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: 1000.0 * (temperature_k - 200.0), temperature_k,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
        }).collect())?)
}

fn boundary(htc: f64) -> resolved::Boundary {
    resolved::Boundary { ambient_k: 300.0, htc_w_m2_k: htc, reservoir: None, source_w: 0.0 }
}

fn parcel(id: u64, column: usize, layer: usize, bed: bool, reverse: bool) -> resolved::DepositedCell {
    let length = 2e-3;
    let width = 4e-4;
    let height = if bed { 1e-3 } else { 2e-4 };
    let left = column as f64 * length;
    let z = if bed { -0.5 * height } else { (layer as f64 + 0.5) * height };
    let endpoints = if reverse { [left+length, left] } else { [left, left+length] };
    let temperature = if bed { 330.0 } else { 500.0 };
    let volume = length * width * height;
    resolved::DepositedCell {
        id: resolved::CellId(id), body_id: if bed { 0 } else { 1 }, segment_id: id,
        material_id: resolved::MaterialId(if bed { 2 } else { 1 }),
        mass_kg: 1000.0 * volume, specific_h_j_kg: 1000.0 * (temperature - 200.0),
        geometry: resolved::CellGeometry { start_m: [endpoints[0], width/2.0, z],
            end_m: [endpoints[1], width/2.0, z], volume_m3: volume,
            surface_area_m2: 2.0 * (length*width + length*height + width*height) },
        boundary: boundary(10.0),
    }
}

fn append_contact(events: &mut Vec<resolved::TimedEvent>, next_id: &mut u64,
    time_s: f64, a: u64, b: u64, area_m2: f64) {
    events.push(resolved::TimedEvent { time_s, event: resolved::Event::Contact(resolved::ThermalLink {
        id: *next_id, a: resolved::CellId(a), b: resolved::CellId(b), area_m2,
        // Declared thermal coefficient only. No healing/strength is inferred.
        conductance_w_k: 12.5 * area_m2,
    }) });
    *next_id += 1;
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "usage: implicit_deposited_strip [segments layers max_dt_s cooling_s]").into());
    }
    let segments: usize = args.first().map_or(Ok(12), |s| s.parse())?;
    let layers: usize = args.get(1).map_or(Ok(3), |s| s.parse())?;
    let max_dt_s: f64 = args.get(2).map_or(Ok(0.25), |s| s.parse())?;
    let cooling_s: f64 = args.get(3).map_or(Ok(5.0), |s| s.parse())?;
    if !(1..=256).contains(&segments) || !(1..=16).contains(&layers)
        || !max_dt_s.is_finite() || !(0.001..=1.0).contains(&max_dt_s)
        || !cooling_s.is_finite() || !(0.0..=100.0).contains(&cooling_s) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "requires segments 1..256, layers 1..16, dt 0.001..1 s, dwell 0..100 s").into());
    }
    let polymer = curve(0x61)?;
    let bed = curve(0x62)?;
    let cell_count = segments * (layers + 1);
    let mut sim = DepositionSimulation::new(vec![
        (MaterialId(1), DepositionMaterial { curve: &polymer, conductivity_w_m_k: 0.2, emissivity: 0.8 }),
        (MaterialId(2), DepositionMaterial { curve: &bed, conductivity_w_m_k: 10.0, emissivity: 0.8 }),
    ], DepositionBudget { max_cells: cell_count, max_contacts: 3*cell_count,
        max_materials: 2, max_events: 5*cell_count, max_steps: 100_000,
        max_operator_solves: 20_000_000 })?;
    let mut events = Vec::new();
    let mut interface = 1;
    // A finite bed is part of the energy ledger. It is not an infinite thermostat.
    for column in 0..segments {
        events.push(resolved::TimedEvent { time_s: 0.0,
            event: resolved::Event::Birth(parcel(1+column as u64, column, 0, true, false)) });
        if column > 0 {
            append_contact(&mut events, &mut interface, 0.0,
                column as u64, 1+column as u64, 4e-7);
        }
    }
    let mut last_birth_s = 0.0;
    for layer in 0..layers {
        let reverse = layer % 2 == 1;
        for step in 0..segments {
            let column = if reverse { segments-1-step } else { step };
            let id = 1 + segments as u64 + (layer*segments + column) as u64;
            last_birth_s = 0.25 * (1 + layer*segments + step) as f64;
            events.push(resolved::TimedEvent { time_s: last_birth_s,
                event: resolved::Event::Birth(parcel(id, column, layer, false, reverse)) });
            let below = if layer == 0 { 1 + column as u64 } else { id - segments as u64 };
            append_contact(&mut events, &mut interface, last_birth_s, below, id, 8e-7);
            if step > 0 {
                let previous = if reverse { id+1 } else { id-1 };
                append_contact(&mut events, &mut interface, last_birth_s, previous, id, 8e-8);
            }
        }
    }
    // An explicit fan change at the final birth, not an inferred controller.
    for id in (1+segments)..=cell_count {
        events.push(resolved::TimedEvent { time_s: last_birth_s,
            event: resolved::Event::Boundary(resolved::CellId(id as u64), boundary(25.0)) });
    }
    let target = last_birth_s + cooling_s;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let started = Instant::now();
    let report = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0x4646_4654, kernel_id: 51,
            tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        sim.run_resolved_until(&cx, &events, EventId(1), target, 300.0,
            DepositionStepConfig { maximum_step_s: max_dt_s,
                solver_energy_tolerance_j: 1e-12, energy_tolerance_j: 1e-7 })
    })?;
    let elapsed_s = started.elapsed().as_secs_f64();
    println!("{{\"model\":\"synthetic-isothermal-implicit-strip\",\"cells\":{},\"contacts\":{},\"events\":{},\"steps\":{},\"operator_solves\":{},\"native_elapsed_s\":{:.6},\"time_s\":{:.6},\"mass_kg\":{:.12e},\"incoming_enthalpy_j\":{:.12e},\"boundary_heat_j\":{:.12e},\"stored_enthalpy_j\":{:.12e},\"residual_j\":{:.12e},\"maximum_biot\":{:.6}}}",
        sim.cells().len(), sim.contacts().len(), report.events, report.steps, report.operator_solves,
        elapsed_s, sim.time_s(), sim.total_mass_kg(), report.incoming_enthalpy_j,
        report.boundary_heat_j, sim.stored_enthalpy_j(), report.energy_residual_j, report.maximum_biot);
    println!("cell,body,segment,birth_s,x_m,y_m,z_m,mass_kg,temperature_K,specific_enthalpy_J_kg");
    for (id, cell) in sim.cells() {
        let p = cell.birth.geometry.position_m;
        println!("{},{},{},{:.6},{:.9},{:.9},{:.9},{:.12e},{:.9},{:.9}", id.0,
            cell.birth.body.0, cell.birth.segment.0, cell.birth_time_s, p[0], p[1], p[2],
            cell.birth.mass_kg, cell.phase.temperature_k(), cell.phase.specific_enthalpy_j_kg());
    }
    Ok(())
}
