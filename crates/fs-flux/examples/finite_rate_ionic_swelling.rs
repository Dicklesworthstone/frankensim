//! Finite-rate ionic and neutral exchange change osmotic swelling trajectories.
//!
//! Run: cargo run -p fs-flux --features continuum --example finite_rate_ionic_swelling
//! All coefficients are illustrative, uncalibrated inputs. Two finite baths
//! share the original virial free energy and original Simo--Taylor mechanics.
//! The second experiment resolves separate sodium, potassium and chloride
//! channels, a neutral permeant, and a retained neutral tracer. Species names
//! describe roles; these coefficients are not fitted electrolyte data.
//! The initial nonequilibrium ion partition is retained. Selective water pores
//! and independent ion channels share the membrane area; no solvent drag,
//! double-layer capacitance, active pump or heat production is inferred.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticPressureVolume,
    electrodiffusion::{IonicAdvanceControls, IonicOsmoticPair},
};
use fs_material::{
    hyper::{
        Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric,
        discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings},
    },
    osmosis::{
        AqueousSolution, HydraulicPermeance, SolutePermeance, donnan::DonnanPotential,
        electrodiffusion::ElectrodiffusiveMembrane, nonideal::VirialSolutionModel,
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = R * T;
const VOLUME: f64 = 1e-6;
const MOLES: f64 = 1e-4;
const AREA: f64 = 1e-4;
const TIME_SCALE: f64 = 1000.0;
const CONCENTRATION: f64 = MOLES / VOLUME;
const B: f64 = 0.8 / CONCENTRATION;

fn mechanics() -> Result<OsmoticPressureVolume, Box<dyn std::error::Error>> {
    let solid = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 1e4,
                alpha: 2.0,
            }],
            volumetric: Volumetric::SimoTaylor {
                kappa: RT * CONCENTRATION,
            },
        },
        4.0,
    )?;
    let bulk = DiscreteBulkPoint::try_new(
        &solid,
        1.0,
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.2,
            maximum_volume_ratio: 2.0,
            maximum_iterations: 96,
            absolute_pressure_tolerance_pa: 1e-8,
            relative_pressure_tolerance: 1e-12,
            absolute_work_tolerance_pa: 1e-8,
            relative_work_tolerance: 1e-12,
        },
    )?;
    Ok(OsmoticPressureVolume::from_bulk_material(
        VOLUME, 0.0, 0.0, bulk,
    )?)
}

fn pair(permeance_scale: f64) -> Result<IonicOsmoticPair, Box<dyn std::error::Error>> {
    let model = Arc::new(VirialSolutionModel::new(
        T,
        18e-6,
        1.0,
        vec![1000.0; 2],
        vec![B, 0.0, 0.0, B],
    )?);
    let mechanics = mechanics()?;
    let cell = |amounts| {
        OsmoticCompartment::new(
            AqueousSolution::nonideal(Arc::clone(&model), VOLUME, amounts)?,
            mechanics,
            0.25 * VOLUME,
        )
    };
    let permeance = permeance_scale * VOLUME / (AREA * TIME_SCALE);
    Ok(IonicOsmoticPair::new(
        [
            cell(vec![1.2 * MOLES, 0.2 * MOLES])?,
            cell(vec![0.8 * MOLES; 2])?,
        ],
        ElectrodiffusiveMembrane::new(
            HydraulicPermeance::new(VOLUME / (AREA * TIME_SCALE * RT * CONCENTRATION))?,
            vec![1, -1],
            vec![
                SolutePermeance::new(0.7 * permeance)?,
                SolutePermeance::new(0.2 * permeance)?,
            ],
        )?,
        [-MOLES, 0.0],
        [0.1 * MOLES, 0.0],
        AREA,
        1e-16,
        1 << 20,
    )?)
}

fn mixture() -> Result<IonicOsmoticPair, Box<dyn std::error::Error>> {
    // Ordered Na+, K+, Cl-, neutral permeant, blocked neutral tracer. The
    // explicit symmetric cross terms belong to the same convex free energy
    // used by water activity, ionic drift, and the total-potential check.
    let interactions = [
        [0.8, 0.15, 0.10, 0.0, 0.0],
        [0.15, 0.5, 0.12, 0.0, 0.0],
        [0.10, 0.12, 0.6, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.4, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.25],
    ];
    let model = Arc::new(VirialSolutionModel::new(
        T,
        18e-6,
        1.0,
        vec![1000.0; 5],
        interactions
            .into_iter()
            .flatten()
            .map(|b| b / CONCENTRATION)
            .collect(),
    )?);
    let mechanics = mechanics()?;
    let cell = |amounts: [f64; 5]| {
        OsmoticCompartment::new(
            AqueousSolution::nonideal(
                Arc::clone(&model),
                VOLUME,
                amounts.into_iter().map(|n| n * MOLES).collect(),
            )?,
            mechanics,
            0.25 * VOLUME,
        )
    };
    let permeance = VOLUME / (AREA * TIME_SCALE);
    Ok(IonicOsmoticPair::new(
        [
            cell([1.0, 0.2, 0.2, 0.25, 0.03])?,
            cell([0.2, 0.6, 0.8, 0.05, 0.01])?,
        ],
        ElectrodiffusiveMembrane::new(
            HydraulicPermeance::new(VOLUME / (AREA * TIME_SCALE * RT * CONCENTRATION))?,
            vec![1, 1, -1, 0, 0],
            [0.3, 0.8, 0.2, 0.25, 0.0]
                .into_iter()
                .map(|p| SolutePermeance::new(p * permeance))
                .collect::<Result<Vec<_>, _>>()?,
        )?,
        [-MOLES, 0.0],
        [0.1 * MOLES, 0.0],
        AREA,
        1e-16,
        1 << 20,
    )?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Finite-rate charged swelling; illustrative uncalibrated coefficients");
    println!("All runs start from the SAME actual ion amounts and solvent volumes.");
    println!(
        "permeance_scale time_s left_mL left_anion_umol voltage_mV delta_potential_J accepted rejected"
    );
    for scale in [0.0, 0.05, 1.0, 20.0] {
        let mut pair = pair(scale)?;
        let mut workspace = pair.adaptive_workspace(1 << 20)?;
        let controls = IonicAdvanceControls {
            maximum_step_s: 100.0,
            relative_tolerance: 2e-5,
            absolute_water_tolerance_m3: 1e-16,
            absolute_amount_tolerance_mol: 1e-13,
            ..Default::default()
        };
        let mut potential_change = 0.0;
        for end in [50.0, 200.0, 1000.0, 5000.0, 20000.0] {
            let report = pair.advance(end - pair.time_s(), controls, &mut workspace, || false)?;
            potential_change += report.total_potential_change_j;
            let voltage = match pair.endpoint().potential {
                DonnanPotential::Finite {
                    voltage_left_minus_right_v,
                    ..
                } => {
                    format!("{:.7}", 1000.0 * voltage_left_minus_right_v)
                }
                DonnanPotential::Indeterminate => "unconstrained".into(),
                DonnanPotential::NegativeInfinite => "-infinity".into(),
                DonnanPotential::PositiveInfinite => "+infinity".into(),
            };
            let left = pair.compartments()[0].solution();
            println!(
                "{scale:.2} {end:.1} {:.9} {:.9} {voltage} {potential_change:.9e} {} {}",
                left.solvent_volume_m3() * 1e6,
                left.solute_moles()[1] * 1e6,
                report.accepted_steps,
                report.rejected_steps,
            );
        }
    }
    println!("Independent Na/K/Cl and neutral exchange; illustrative uncalibrated mixture");
    println!(
        "time_s left_mL Na_umol K_umol Cl_umol neutral_umol retained_umol delta_potential_J accepted rejected"
    );
    let mut mixed = mixture()?;
    let mut workspace = mixed.adaptive_workspace(1 << 20)?;
    let controls = IonicAdvanceControls {
        maximum_step_s: 100.0,
        relative_tolerance: 2e-5,
        absolute_water_tolerance_m3: 1e-16,
        absolute_amount_tolerance_mol: 1e-13,
        ..Default::default()
    };
    let mut potential_change = 0.0;
    for end in [50.0, 200.0, 1000.0, 5000.0, 20000.0] {
        let report = mixed.advance(end - mixed.time_s(), controls, &mut workspace, || false)?;
        potential_change += report.total_potential_change_j;
        let left = mixed.compartments()[0].solution();
        let amounts = left.solute_moles();
        println!(
            "{end:.1} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {potential_change:.9e} {} {}",
            left.solvent_volume_m3() * 1e6,
            amounts[0] * 1e6,
            amounts[1] * 1e6,
            amounts[2] * 1e6,
            amounts[3] * 1e6,
            amounts[4] * 1e6,
            report.accepted_steps,
            report.rejected_steps,
        );
    }
    Ok(())
}
