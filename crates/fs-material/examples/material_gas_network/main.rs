//! Four source-backed volumes with two inlets and an outlet at a mixing node.
//! Reuse the original gas-source loader and the shared snapshot nozzle solver.
#[path = "../material_gas_spring/source.rs"]
mod source;

use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::{ClosedGasParcel,
    transfer::{snapshot::GasSnapshotSolve, nozzle::{GasNozzle, network::snapshot::{
        SnapshotNozzleLink, SnapshotNozzleStep, step_candidate,
    }}}}};
use std::{collections::BTreeMap, error::Error, io::Read};
type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasSnapshotSolve {
    GasSnapshotSolve { maximum_withdrawn_fraction: 0.05,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }
}
fn read_source(path: &str) -> Result<String> {
    let mut text = String::new();
    std::fs::File::open(path)?.take(65537).read_to_string(&mut text)?;
    if text.len() > 65536 { return Err("synthetic mixture input exceeds 64 KiB".into()); }
    Ok(text)
}
fn initial(a: &str, b: &str, provenance: &str) -> Result<BTreeMap<u64, ClosedGasParcel>> {
    let a = source::load(a, provenance)?;
    let b = source::load(b, provenance)?;
    let mut parcels = BTreeMap::new();
    for (id, model, t, p) in [(0, a, 500.0, 300000.0), (1, b.clone(), 450.0, 220000.0),
        (2, b.clone(), 400.0, 140000.0), (3, b, 350.0, 100000.0)] {
        parcels.insert(id, ClosedGasParcel::new(model, source::point(t, p)?, 0.001)?);
    }
    Ok(parcels)
}
fn links(area_m2: f64) -> Vec<SnapshotNozzleLink> {
    [(10, 0, 2), (20, 1, 2), (30, 2, 3)].into_iter().map(|(id, a, b)| SnapshotNozzleLink {
        id, parcels: [a, b], nozzle: GasNozzle { area_m2, discharge_coefficient: 0.8 }, nozzle_solve: energy(),
    }).collect()
}
fn advance(parcels: &BTreeMap<u64, ClosedGasParcel>, area: f64, dt: f64)
    -> Result<(BTreeMap<u64, ClosedGasParcel>, SnapshotNozzleStep)> {
    let endpoints = parcels.keys().map(|&id| (id, energy())).collect();
    Ok(step_candidate(parcels, &links(area), dt, &endpoints, controls(), || false)?.into_parts())
}
fn energy_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    next.mass_kg()*(next.state().specific_internal_energy_j_kg()-old.state().specific_internal_energy_j_kg())
        +(next.mass_kg()-old.mass_kg())*old.state().specific_internal_energy_j_kg()
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 2 || args.len() > 5 {
        return Err("usage: material_gas_network A.tsv B.tsv [steps=200] [dt_s=0.0002] [area_m2=0.000001]".into());
    }
    let steps: usize = args.get(2).map_or(Ok(200), |s| s.parse())?;
    let dt: f64 = args.get(3).map_or(Ok(0.0002), |s| s.parse())?;
    let area: f64 = args.get(4).map_or(Ok(1e-6), |s| s.parse())?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt <= 0.0
        || !area.is_finite() || area < 0.0 || !(dt*steps as f64).is_finite() {
        return Err("bounded step count, positive finite dt and nonnegative finite area required".into());
    }
    let provenance = format!("one synthetic source set: {} + {}", args[0], args[1]);
    let original = initial(&read_source(&args[0])?, &read_source(&args[1])?, &provenance)?;
    let mut parcels = original.clone();
    println!("time_s,parcel_id,mass_kg,temperature_k,pressure_pa,sound_speed_m_s,impedance_pa_s_m,viscosity_pa_s,conductivity_w_m_k,incoming_mass_kg,outgoing_mass_kg,network_energy_error_j,step_species_error_kg,choked_links");
    for i in 0..=steps {
        let report = if i == 0 { None } else {
            let (next, report) = advance(&parcels, area, dt)?; parcels = next; Some(report)
        };
        let drift: f64 = parcels.iter().map(|(id, p)| energy_change(&original[id], p)).sum();
        for (id, parcel) in &parcels {
            let gas = parcel.state().state();
            let account = report.as_ref().and_then(|r| r.exchange.parcels.iter().find(|p| &p.id == id));
            println!("{:.9},{},{:.12e},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{}",
                i as f64*dt, id, parcel.mass_kg(), gas.temperature, gas.pressure, gas.sound_speed,
                gas.characteristic_impedance, gas.dynamic_viscosity, gas.thermal_conductivity,
                account.map_or(0.0, |a| a.incoming_mass_kg), account.map_or(0.0, |a| a.outgoing_mass_kg),
                drift, report.as_ref().map_or(0.0, |r| r.exchange.maximum_species_mass_residual_kg),
                report.as_ref().map_or(0, |r| r.links.iter().filter(|l| l.flow.choked).count()));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const A: &str = include_str!("../material_gas_spring/mixture-a.tsv");
    const B: &str = include_str!("../material_gas_spring/mixture-b.tsv");
    #[test]
    fn actual_source_files_drive_branching_conservative_histories() {
        let mut outcomes = Vec::new();
        for (a,b) in [(A,B),(B,A)] {
            let original = initial(a,b,"shared synthetic test sources").unwrap();
            let mut parcels = original.clone();
            for _ in 0..200 {
                let (next, report) = advance(&parcels, 1e-6, 0.0002).unwrap();
                assert!(report.exchange.total_energy_residual_j.abs() < 1e-8);
                parcels = next;
            }
            let drift: f64 = parcels.iter().map(|(id,p)| energy_change(&original[id],p)).sum();
            assert!(drift.abs() < 1e-6);
            assert_ne!(parcels[&2].state().state().specific_gas_constant, original[&2].state().state().specific_gas_constant);
            outcomes.push(parcels[&3].state().state().sound_speed);
            let (closed, _) = advance(&parcels, 0.0, 0.0002).unwrap();
            assert_eq!(closed, parcels);
        }
        assert_ne!(outcomes[0], outcomes[1]);
    }
}
