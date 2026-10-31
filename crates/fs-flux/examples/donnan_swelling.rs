//! Finite charged-compartment swelling in the ideal, fast-ion limit.
//!
//! Run with `cargo run -p fs-flux --example donnan_swelling`.
//! These explicit coefficients are an illustrative model, not measured material
//! data. The left side contains 1 micromole of trapped particles carrying
//! 100 micromoles of negative charge equivalents. That charge is not counted
//! again as 100 micromoles of osmotic particles. Both salt baths are finite.
//! Linear chamber compliance arrests swelling; no infinite salt reservoir,
//! arbitrary ionic reflection factor or imposed electric potential is used.

use fs_flux::osmosis::donnan::{DonnanOsmoticPair, DonnanStepControls};
use fs_flux::osmosis::{OsmoticCompartment, OsmoticPressureVolume};
use fs_material::osmosis::donnan::{DonnanPotential, IdealDonnan};
use fs_material::osmosis::{HydraulicPermeance, IdealDiluteSolution};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const TEMPERATURE: f64 = 298.15;
const INITIAL_VOLUME: f64 = 1e-6;
const SALT_UNIT: f64 = 1e-4;
const TRAPPED: f64 = 1e-6;
const ELASTANCE: f64 = 2.5e11;

fn compartment(amounts: Vec<f64>) -> Result<OsmoticCompartment, Box<dyn std::error::Error>> {
    Ok(OsmoticCompartment::new(
        IdealDiluteSolution::new(TEMPERATURE, 18e-6, INITIAL_VOLUME, amounts)?,
        OsmoticPressureVolume::new(INITIAL_VOLUME, 0.0, ELASTANCE)?,
        0.05 * INITIAL_VOLUME,
    )?)
}

// Independent monovalent algebra: n_left_minus=a*N, n_left_plus=(a+1)*N,
// n_right_minus=n_right_plus=(1-a)*N. Equal electrochemical potentials imply
// a*(a+1)/VL²=(1-a)²/VR². Use its stable positive quadratic root.
fn equilibrium_left_volume() -> f64 {
    let driving_pressure = |left: f64| {
        let right = 2.0 * INITIAL_VOLUME - left;
        let l = left / INITIAL_VOLUME;
        let r = right / INITIAL_VOLUME;
        let b = r * r + 2.0 * l * l;
        let a = 2.0 * l * l / (b + (b * b + 4.0 * (r * r - l * l) * l * l).sqrt());
        let osmotic = R
            * TEMPERATURE
            * ((SALT_UNIT * (2.0 * a + 1.0) + TRAPPED) / left
                - 2.0 * SALT_UNIT * (1.0 - a) / right);
        2.0 * ELASTANCE * (left - INITIAL_VOLUME) - osmotic
    };
    let (mut low, mut high) = (INITIAL_VOLUME, 1.95 * INITIAL_VOLUME);
    assert!(driving_pressure(low) < 0.0 && driving_pressure(high) > 0.0);
    for _ in 0..80 {
        let midpoint = low + 0.5 * (high - low);
        if driving_pressure(midpoint) < 0.0 {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    low + 0.5 * (high - low)
}

fn free_energy(pair: &DonnanOsmoticPair) -> f64 {
    let mut mixing = 0.0;
    let mut mechanical = 0.0;
    for (side, cell) in pair.compartments().iter().enumerate() {
        let volume = cell.solution().solvent_volume_m3();
        for &amount in cell
            .solution()
            .solute_moles()
            .iter()
            .chain(std::iter::once(&pair.trapped_particles_mol()[side]))
        {
            if amount > 0.0 {
                // c_ref=1 mol/m³; its fixed-inventory constant cancels in differences.
                mixing += amount * ((amount / volume).ln() - 1.0);
            }
        }
        mechanical += 0.5 * ELASTANCE * (volume - INITIAL_VOLUME).powi(2);
    }
    R * TEMPERATURE * mixing + mechanical
}

fn print_state(pair: &DonnanOsmoticPair) {
    let state = pair.equilibrium();
    let voltage = match state.potential {
        DonnanPotential::Finite {
            voltage_left_minus_right_v,
            ..
        } => voltage_left_minus_right_v,
        _ => panic!("the supplied two-ion finite-salt fixture has a finite potential"),
    };
    let left = &pair.compartments()[0];
    let right = &pair.compartments()[1];
    println!(
        "{:8.1} {:12.6} {:12.6} {:12.6} {:14.6}",
        pair.time_s(),
        left.solution().solvent_volume_m3() * 1e6,
        right.solution().solvent_volume_m3() * 1e6,
        voltage * 1e3,
        state.osmotic_pressure_pa[0] - state.osmotic_pressure_pa[1],
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut pair = DonnanOsmoticPair::new(
        [
            compartment(vec![1.5 * SALT_UNIT, 0.5 * SALT_UNIT])?,
            compartment(vec![0.5 * SALT_UNIT, 0.5 * SALT_UNIT])?,
        ],
        IdealDonnan::new(vec![1, -1])?,
        [-SALT_UNIT, 0.0],
        [TRAPPED, 0.0],
        1e-4,
        HydraulicPermeance::new(1e-11)?,
        1 << 20,
    )?;
    let initial_energy = free_energy(&pair);
    let mut pressure_work = 0.0;
    let mut maximum_water_residual = 0.0_f64;
    let mut maximum_ion_residual = 0.0_f64;
    let mut maximum_charge_residual = 0.0_f64;
    let mut maximum_energy_increase = 0.0_f64;
    let mut previous_energy = initial_energy;
    println!("Ideal fast-ion Donnan swelling; illustrative uncalibrated coefficients");
    println!("  time_s    left_mL     right_mL    potential_mV  osmotic_drop_Pa");
    print_state(&pair);
    for end in [100.0, 1000.0, 5000.0, 12000.0, 24000.0] {
        while pair.time_s() < end {
            let report = pair.step(
                (end - pair.time_s()).min(5.0),
                DonnanStepControls::default(),
                || false,
            )?;
            pressure_work += report.pressure_work_into_solution_j;
            maximum_water_residual =
                maximum_water_residual.max(report.water_balance_residual_m3.abs());
            for &residual in &report.solute_balance_residual_mol[..report.species_count] {
                maximum_ion_residual = maximum_ion_residual.max(residual.abs());
            }
            for &residual in &report.equilibrium.charge_residual_mol {
                maximum_charge_residual = maximum_charge_residual.max(residual.abs());
            }
            let energy = free_energy(&pair);
            maximum_energy_increase = maximum_energy_increase.max(energy - previous_energy);
            previous_energy = energy;
        }
        print_state(&pair);
    }
    let left = pair.compartments()[0].solution().solvent_volume_m3();
    let exact = equilibrium_left_volume();
    println!("Independent equilibrium left volume: {:.9} mL", exact * 1e6);
    println!("Remaining volume error: {:.3e} m³", left - exact);
    println!(
        "Mixing plus chamber energy change: {:.9e} J",
        free_energy(&pair) - initial_energy
    );
    println!("Largest per-step energy increase: {maximum_energy_increase:.3e} J");
    println!("Mechanical work into solutions: {pressure_work:.9e} J");
    println!(
        "Maximum step residuals: water={maximum_water_residual:.3e} m³, ion={maximum_ion_residual:.3e} mol, charge={maximum_charge_residual:.3e} mol equivalents"
    );
    Ok(())
}
