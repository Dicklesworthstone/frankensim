//! Finite nonideal charged-bath swelling against original nonlinear mechanics.
//!
//! Run: cargo run -p fs-flux --features continuum --example nonideal_donnan_swelling
//! The coefficients below are authored illustrative inputs, not a calibration
//! of a particular salt or gel. Both baths are finite. Fixed charge, separate
//! trapped particles, retained dry volume and mobile ions remain distinct.
//! Ion activities and solvent pressure share one original virial free energy;
//! the Simo--Taylor chamber law supplies its own pressure and exact storage.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticPressureVolume,
    donnan::{DonnanAdvanceControls, DonnanOsmoticPair},
};
use fs_material::{
    hyper::{
        Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric,
        discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings},
    },
    osmosis::{
        AqueousSolution, HydraulicPermeance,
        donnan::{DonnanPotential, VirialDonnan},
        nonideal::VirialSolutionModel,
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = R * T;
const V: f64 = 1e-6;
const DRY_VOLUME: f64 = 0.2e-6;
const MATERIAL_VOLUME: f64 = V + DRY_VOLUME;
const SALT: f64 = 1e-4;
const TRAPPED: f64 = 1e-5;
const B: f64 = 0.008;
const BULK_MODULUS: f64 = 2.5e5;

// Independent scalar ion oracle: charge gives nL=SALT*[1+a,a] and
// nR=SALT*[1-a,1-a]. Summing the two activity equations eliminates voltage.
fn partition(left: f64) -> f64 {
    let right = 2.0 * V - left;
    let residual = |a: f64| {
        ((1.0 + a) * SALT / left).ln() + (a * SALT / left).ln()
            - 2.0 * ((1.0 - a) * SALT / right).ln()
            + B * SALT * ((1.0 + 2.0 * a) / left - 2.0 * (1.0 - a) / right)
    };
    bisect(1e-12, 1.0 - 1e-12, residual)
}
fn bisect(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lo) < 0.0 && f(hi) > 0.0);
    for _ in 0..80 {
        let mid = lo + 0.5 * (hi - lo);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}
fn chamber_pressure(volume: f64) -> f64 {
    let j = (volume + DRY_VOLUME) / MATERIAL_VOLUME;
    0.5 * BULK_MODULUS * (j - 1.0 / j)
}
fn osmotic(volume: f64, amounts: [f64; 2], trapped: f64) -> f64 {
    RT * (amounts
        .into_iter()
        .map(|n| {
            let c = n / volume;
            c + 0.5 * B * c * c
        })
        .sum::<f64>()
        + trapped / volume)
}
fn equilibrium_volume() -> f64 {
    bisect(V, 1.45 * V, |left| {
        let right = 2.0 * V - left;
        let a = partition(left);
        chamber_pressure(left)
            - chamber_pressure(right)
            - osmotic(left, [SALT * (1.0 + a), SALT * a], TRAPPED)
            + osmotic(right, [SALT * (1.0 - a); 2], 0.0)
    })
}
fn energy(pair: &DonnanOsmoticPair) -> f64 {
    pair.compartments()
        .iter()
        .enumerate()
        .map(|(side, cell)| {
            let volume = cell.solution().solvent_volume_m3();
            let mixing = cell
                .solution()
                .solute_moles()
                .iter()
                .map(|&n| n * ((n / volume).ln() - 1.0) + 0.5 * B * n * n / volume)
                .sum::<f64>();
            let trapped = if side == 0 {
                TRAPPED * ((TRAPPED / volume).ln() - 1.0)
            } else {
                0.0
            };
            let j = (volume + DRY_VOLUME) / MATERIAL_VOLUME;
            RT * (mixing + trapped)
                + 0.25 * MATERIAL_VOLUME * BULK_MODULUS * (j * j - 1.0 - 2.0 * j.ln())
        })
        .sum()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let solution = Arc::new(VirialSolutionModel::new(
        T,
        18e-6,
        1.0,
        vec![1000.0; 2],
        vec![B, 0.0, 0.0, B],
    )?);
    let material = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 1e4,
                alpha: 2.0,
            }],
            volumetric: Volumetric::SimoTaylor {
                kappa: BULK_MODULUS,
            },
        },
        4.0,
    )?;
    let point = DiscreteBulkPoint::try_new(
        &material,
        1.0,
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.4,
            maximum_volume_ratio: 2.5,
            maximum_iterations: 96,
            absolute_pressure_tolerance_pa: 1e-8,
            relative_pressure_tolerance: 1e-12,
            absolute_work_tolerance_pa: 1e-8,
            relative_work_tolerance: 1e-12,
        },
    )?;
    let pressure = OsmoticPressureVolume::from_bulk_material(V, DRY_VOLUME, 0.0, point)?;
    let cell = |amounts| {
        OsmoticCompartment::new(
            AqueousSolution::nonideal(Arc::clone(&solution), V, amounts)?,
            pressure,
            0.05 * V,
        )
    };
    let mut pair = DonnanOsmoticPair::new_nonideal(
        [
            cell(vec![1.5 * SALT, 0.5 * SALT])?,
            cell(vec![0.5 * SALT; 2])?,
        ],
        VirialDonnan::new(vec![1, -1], Arc::clone(&solution), 2e-12)?,
        [-SALT, 0.0],
        [TRAPPED, 0.0],
        1e-4,
        HydraulicPermeance::new(1e-11)?,
        1 << 20,
    )?;
    let mut workspace = pair.adaptive_workspace(1 << 20)?;
    let controls = DonnanAdvanceControls {
        maximum_step_s: 500.0,
        relative_tolerance: 2e-5,
        absolute_water_tolerance_m3: 1e-16,
        absolute_amount_tolerance_mol: 1e-13,
        ..Default::default()
    };
    let initial_energy = energy(&pair);
    let mut work = 0.0;
    let mut previous_energy = initial_energy;
    println!("Nonideal fast-ion swelling; illustrative uncalibrated coefficients");
    println!(
        "time_s left_mL right_mL voltage_mV osmotic_drop_Pa delta_potential_J accepted rejected"
    );
    for end in [100.0, 500.0, 2000.0, 6000.0, 12000.0, 25000.0] {
        let report = pair.advance(end - pair.time_s(), controls, &mut workspace, || false)?;
        let state = pair.equilibrium();
        let voltage = match state.potential {
            DonnanPotential::Finite {
                voltage_left_minus_right_v,
                ..
            } => voltage_left_minus_right_v,
            _ => return Err("the finite two-ion example requires a finite voltage".into()),
        };
        let current_energy = energy(&pair);
        assert!(current_energy <= previous_energy + 1e-10);
        assert!(report.water_balance_residual_m3.abs() <= 1e-15);
        assert!(
            report.solute_balance_residual_mol[..2]
                .iter()
                .all(|r| r.abs() <= 1e-13)
        );
        work += report.pressure_work_into_solution_j;
        previous_energy = current_energy;
        println!(
            "{end:.1} {:.9} {:.9} {:.6} {:.6} {:.9e} {} {}",
            pair.compartments()[0].solution().solvent_volume_m3() * 1e6,
            pair.compartments()[1].solution().solvent_volume_m3() * 1e6,
            voltage * 1e3,
            state.osmotic_pressure_pa[0] - state.osmotic_pressure_pa[1],
            current_energy - initial_energy,
            report.accepted_steps,
            report.rejected_steps
        );
    }
    let expected = equilibrium_volume();
    println!(
        "Independent equilibrium left volume: {:.9} mL",
        expected * 1e6
    );
    println!(
        "Remaining volume error: {:.6e} m3",
        pair.compartments()[0].solution().solvent_volume_m3() - expected
    );
    println!("Integrated mechanical work into solutions: {work:.9e} J");
    Ok(())
}
