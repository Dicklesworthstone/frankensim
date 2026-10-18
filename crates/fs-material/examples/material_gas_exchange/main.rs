//! Two finite material gas volumes exchanging species and total enthalpy.
//! Synthetic source inputs use the existing gas-spring loader, not another schema.
#[path = "../material_gas_spring/source.rs"]
mod source;

use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::{ClosedGasParcel,
    transfer::{GasMassTransferSolve, nozzle::{GasNozzle, GasNozzleStepSolve}}}};
use std::{error::Error, io::Read};
type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn controls() -> GasNozzleStepSolve {
    let energy = GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 };
    GasNozzleStepSolve { nozzle: energy, maximum_donor_fraction: 0.05,
        transfer: GasMassTransferSolve { donor: energy, receiver: energy,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn read_source(path: &str) -> Result<String> {
    let mut text = String::new();
    std::fs::File::open(path)?.take(65537).read_to_string(&mut text)?;
    if text.len() > 65536 { return Err("synthetic mixture input exceeds 64 KiB".into()); }
    Ok(text)
}
fn initial(left: &str, right: &str, origin: &str) -> Result<(ClosedGasParcel, ClosedGasParcel)> {
    // Both synthetic files explicitly form ONE source set. Preserve both paths
    // in their common provenance; inconsistent same-species laws still refuse
    // during transfer. Production source-pack joins remain a separate API.
    let a = source::load(left, origin)?;
    let b = source::load(right, origin)?;
    Ok((ClosedGasParcel::new(a, source::point(500.0, 300000.0)?, 0.001)?,
        ClosedGasParcel::new(b, source::point(400.0, 100000.0)?, 0.001)?))
}

fn energy_change(initial: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    next.mass_kg() * (next.state().specific_internal_energy_j_kg()
        - initial.state().specific_internal_energy_j_kg())
        + (next.mass_kg() - initial.mass_kg()) * initial.state().specific_internal_energy_j_kg()
}

fn advance(
    a: &ClosedGasParcel, b: &ClosedGasParcel, area_m2: f64, dt: f64,
) -> Result<(ClosedGasParcel, ClosedGasParcel, f64, bool)> {
    let nozzle = GasNozzle { area_m2, discharge_coefficient: 0.8 };
    if a.state().state().pressure >= b.state().state().pressure {
        let (trial, flow) = a.through_nozzle_to(b, nozzle, dt, controls(), || false)?;
        let (a, b, _) = trial.into_parts();
        Ok((a, b, flow.mass_flow_kg_s, flow.choked))
    } else {
        let (trial, flow) = b.through_nozzle_to(a, nozzle, dt, controls(), || false)?;
        let (b, a, _) = trial.into_parts();
        Ok((a, b, -flow.mass_flow_kg_s, flow.choked))
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 2 || args.len() > 5 {
        return Err("usage: material_gas_exchange LEFT.tsv RIGHT.tsv [steps=200] [dt_s=0.0002] [area_m2=0.000001]".into());
    }
    let steps: usize = args.get(2).map_or(Ok(200), |s| s.parse())?;
    let dt: f64 = args.get(3).map_or(Ok(0.0002), |s| s.parse())?;
    let area: f64 = args.get(4).map_or(Ok(1e-6), |s| s.parse())?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt <= 0.0
        || !area.is_finite() || area < 0.0 || !(dt * steps as f64).is_finite() {
        return Err("bounded step count, positive finite dt and nonnegative finite area required".into());
    }
    let provenance = format!("paired synthetic source set: {} + {}", args[0], args[1]);
    let (mut a, mut b) = initial(&read_source(&args[0])?, &read_source(&args[1])?, &provenance)?;
    let (initial_a, initial_b) = (a.clone(), b.clone());
    let mut transferred = 0.0;
    println!("time_s,left_mass_kg,right_mass_kg,left_temperature_k,right_temperature_k,left_pressure_pa,right_pressure_pa,left_sound_speed_m_s,right_sound_speed_m_s,right_specific_gas_constant_j_kg_k,right_impedance_pa_s_m,right_viscosity_pa_s,right_conductivity_w_m_k,mass_rate_left_to_right_kg_s,choked,accumulated_mass_left_to_right_kg,pair_energy_error_j,pair_mass_error_kg");
    for i in 0..=steps {
        let (rate, choked) = if i == 0 { (0.0, false) } else {
            let (next_a, next_b, rate, choked) = advance(&a, &b, area, dt)?;
            a = next_a; b = next_b; transferred += rate*dt;
            (rate, choked)
        };
        let left = a.state().state(); let right = b.state().state();
        let energy_error = energy_change(&initial_a, &a) + energy_change(&initial_b, &b);
        let mass_error = (a.mass_kg()-initial_a.mass_kg()) + (b.mass_kg()-initial_b.mass_kg());
        println!("{:.9},{:.12e},{:.12e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e},{:.12e}",
            i as f64*dt, a.mass_kg(), b.mass_kg(), left.temperature, right.temperature,
            left.pressure, right.pressure, left.sound_speed, right.sound_speed,
            right.specific_gas_constant, right.characteristic_impedance,
            right.dynamic_viscosity, right.thermal_conductivity, rate, u8::from(choked),
            transferred, energy_error, mass_error);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const A: &str = include_str!("../material_gas_spring/mixture-a.tsv");
    const B: &str = include_str!("../material_gas_spring/mixture-b.tsv");
    #[test]
    fn actual_data_files_drive_conservative_flow_and_data_only_swaps() {
        let mut final_states = Vec::new();
        for (left, right) in [(A, B), (B, A)] {
            let (mut a, mut b) = initial(left, right, "original paired synthetic fixtures").unwrap();
            let (initial_a, initial_b) = (a.clone(), b.clone());
            for _ in 0..200 {
                let (next_a, next_b, rate, _) = advance(&a, &b, 1e-6, 0.0002).unwrap();
                assert!(rate > 0.0);
                a = next_a; b = next_b;
                assert!((energy_change(&initial_a, &a)+energy_change(&initial_b, &b)).abs() < 1e-7);
            }
            assert_ne!(b.state().state().specific_gas_constant, initial_b.state().state().specific_gas_constant);
            final_states.push((a.state().state().pressure, b.state().state().sound_speed));
        }
        assert!((final_states[0].0-final_states[1].0).abs() > 1000.0);
        assert_ne!(final_states[0].1, final_states[1].1);
    }
    #[test]
    fn closed_path_and_failed_large_step_preserve_both_accepted_volumes() {
        let (a, b) = initial(A, B, "original paired synthetic fixtures").unwrap();
        let saved = (a.clone(), b.clone());
        let (next_a, next_b, rate, _) = advance(&a, &b, 0.0, 0.0002).unwrap();
        assert_eq!((next_a, next_b), saved);
        assert_eq!(rate, 0.0);
        assert!(advance(&a, &b, 1e-6, 1000.0).is_err());
        assert_eq!((a, b), saved);
    }
}
