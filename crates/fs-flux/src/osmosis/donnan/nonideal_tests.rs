use super::*;
use crate::osmosis::{IdealDiluteSolution, OsmoticPressureVolume};
use fs_material::{
    hyper::{
        Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric,
        discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings},
    },
    osmosis::{donnan::DonnanPotential, nonideal::VirialSolutionModel},
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = R * T;

fn model(b: f64) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, 18e-6, 1.0, vec![8.0; 2], vec![b, 0.0, 0.0, b]).unwrap())
}
fn cell(solution: AqueousSolution, gauge: f64) -> OsmoticCompartment {
    let material = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 3.0,
                alpha: 2.0,
            }],
            volumetric: Volumetric::SimoTaylor { kappa: RT },
        },
        4.0,
    )
    .unwrap();
    let point = DiscreteBulkPoint::try_new(
        &material,
        1.0,
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.2,
            maximum_volume_ratio: 2.6,
            maximum_iterations: 96,
            absolute_pressure_tolerance_pa: 1e-10,
            relative_pressure_tolerance: 1e-13,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-13,
        },
    )
    .unwrap();
    OsmoticCompartment::new(
        solution,
        OsmoticPressureVolume::from_bulk_material(1.0, 0.0, gauge, point).unwrap(),
        0.05,
    )
    .unwrap()
}
fn cells(source: Option<&Arc<VirialSolutionModel>>, gauge: f64) -> [OsmoticCompartment; 2] {
    [[1.5, 0.5], [0.5, 0.5]].map(|amounts| {
        let solution = match source {
            Some(source) => {
                AqueousSolution::nonideal(Arc::clone(source), 1.0, amounts.to_vec()).unwrap()
            }
            None => IdealDiluteSolution::new(T, 18e-6, 1.0, amounts.to_vec())
                .unwrap()
                .into(),
        };
        cell(solution, gauge)
    })
}
fn nonideal(b: f64, gauge: f64) -> DonnanOsmoticPair {
    let source = model(b);
    DonnanOsmoticPair::new_nonideal(
        cells(Some(&source), gauge),
        VirialDonnan::new(vec![1, -1], source, 1e-11).unwrap(),
        [-1.0, 0.0],
        [0.1, 0.0],
        1.0,
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        1 << 20,
    )
    .unwrap()
}
fn step_controls() -> DonnanStepControls {
    DonnanStepControls {
        absolute_water_tolerance_m3: 1e-13,
        relative_water_tolerance: 1e-12,
        ..DonnanStepControls::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn root(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
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
// Charge and species constraints leave one ion coordinate:
// nL=[1+t,t], nR=[1-t,1-t]. Summing the two electrochemical equations
// eliminates the membrane potential without calling a production ion solver.
fn ion_root(x: f64, b: f64) -> f64 {
    let (vl, vr) = (1.0 - x, 1.0 + x);
    root(1e-12, 1.0 - 1e-12, |t| {
        (t * (1.0 + t) * vr * vr / ((1.0 - t).powi(2) * vl * vl)).ln()
            + b * ((2.0 * t + 1.0) / vl - 2.0 * (1.0 - t) / vr)
    })
}
fn osmotic(amounts: &[f64], volume: f64, b: f64, trapped: f64) -> f64 {
    trapped / volume
        + amounts
            .iter()
            .map(|n| {
                let c = n / volume;
                c + 0.5 * b * c * c
            })
            .sum::<f64>()
}
fn bulk_energy(volume: f64) -> f64 {
    0.25 * RT * (volume * volume - 1.0 - 2.0 * volume.ln())
}
fn cell_energy(volume: f64, amounts: &[f64], b: f64, trapped: f64) -> f64 {
    let mix = amounts
        .iter()
        .map(|n| n * ((n / volume).ln() - 1.0) + 0.5 * b * n * n / volume)
        .sum::<f64>();
    let retained = if trapped == 0.0 {
        0.0
    } else {
        trapped * ((trapped / volume).ln() - 1.0)
    };
    RT * (mix + retained) + bulk_energy(volume)
}
fn energy(pair: &DonnanOsmoticPair, b: f64) -> f64 {
    pair.compartments()
        .iter()
        .enumerate()
        .map(|(side, cell)| {
            cell_energy(
                cell.solution().solvent_volume_m3(),
                cell.solution().solute_moles(),
                b,
                pair.trapped_particles_mol()[side],
            )
        })
        .sum()
}
fn snapshot(pair: &DonnanOsmoticPair) -> Vec<u64> {
    let mut bits = vec![pair.time_s().to_bits()];
    for c in pair.compartments() {
        bits.push(c.solution().solvent_volume_m3().to_bits());
        bits.extend(c.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    let e = pair.equilibrium();
    bits.extend(e.osmotic_pressure_pa.map(f64::to_bits));
    bits.extend(e.mobile_osmoles_mol.map(f64::to_bits));
    bits.extend(e.charge_residual_mol.map(f64::to_bits));
    bits.push(e.iterations as u64);
    bits.push(
        e.osmotic_pressure_drop_slope_pa_per_m3
            .map_or(u64::MAX, f64::to_bits),
    );
    match e.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            voltage_left_minus_right_v,
        } => bits.extend([
            0,
            dimensionless_left_minus_right.to_bits(),
            voltage_left_minus_right_v.to_bits(),
        ]),
        DonnanPotential::NegativeInfinite => bits.push(1),
        DonnanPotential::PositiveInfinite => bits.push(2),
        DonnanPotential::Indeterminate => bits.push(3),
    }
    bits
}

#[test]
fn zero_interaction_recovers_ideal_water_and_ion_trajectories() {
    let mut ideal = DonnanOsmoticPair::new(
        cells(None, 0.0),
        IdealDonnan::new(vec![1, -1]).unwrap(),
        [-1.0, 0.0],
        [0.1, 0.0],
        1.0,
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        1 << 20,
    )
    .unwrap();
    let mut virial = nonideal(0.0, 0.0);
    for dt in [0.03, 0.12, 0.35] {
        let expected = ideal.step(dt, step_controls(), || false).unwrap();
        let actual = virial.step(dt, step_controls(), || false).unwrap();
        near(
            actual.solvent_transfer_m3,
            expected.solvent_transfer_m3,
            2e-11,
        );
        near(
            actual.pressure_work_into_solution_j,
            expected.pressure_work_into_solution_j,
            2e-7,
        );
        for side in 0..2 {
            let a = virial.compartments()[side].solution();
            let b = ideal.compartments()[side].solution();
            assert!(!a.is_ideal());
            near(a.solvent_volume_m3(), b.solvent_volume_m3(), 3e-11);
            for (&a, &b) in a.solute_moles().iter().zip(b.solute_moles()) {
                near(a, b, 3e-11);
            }
            near(
                actual.equilibrium.osmotic_pressure_pa[side],
                expected.equilibrium.osmotic_pressure_pa[side],
                2e-7,
            );
        }
        near(actual.water_balance_residual_m3, 0.0, 1e-14);
        for residual in &actual.solute_balance_residual_mol[..2] {
            near(*residual, 0.0, 1e-14);
        }
    }
}

#[test]
fn virial_swelling_matches_independent_ion_water_and_energy_oracles_at_any_gauge() {
    let b = 0.8;
    let expected = root(-0.4, 0.0, |x| {
        let (vl, vr, t) = (1.0 - x, 1.0 + x, ion_root(x, b));
        let pressure = 0.5 * (vl - 1.0 / vl - vr + 1.0 / vr);
        let osmotic_drop = osmotic(&[1.0 + t, t], vl, b, 0.1) - osmotic(&[1.0 - t; 2], vr, b, 0.0);
        x - 0.5 * (pressure - osmotic_drop)
    });
    near(expected, -0.13086973250866368, 2e-14);
    let (vl, vr, t, initial_t) = (
        1.0 - expected,
        1.0 + expected,
        ion_root(expected, b),
        ion_root(0.0, b),
    );
    near(t, 0.39142252229687735, 2e-14);
    let exact_work = -bulk_energy(vl) - bulk_energy(vr);
    let expected_energy = cell_energy(vl, &[1.0 + t, t], b, 0.1)
        + cell_energy(vr, &[1.0 - t; 2], b, 0.0)
        - cell_energy(1.0, &[1.0 + initial_t, initial_t], b, 0.1)
        - cell_energy(1.0, &[1.0 - initial_t; 2], b, 0.0);
    near(exact_work, -42.640707594922674, 2e-10);
    near(expected_energy, -148.84253501252635, 2e-9);
    let mut pair = nonideal(b, 0.0);
    let before = energy(&pair, b);
    let report = pair.step(0.5, step_controls(), || false).unwrap();
    near(report.solvent_transfer_m3, expected, 2e-11);
    near(
        pair.compartments()[0].solution().solute_moles()[0],
        1.0 + t,
        2e-11,
    );
    near(
        pair.compartments()[0].solution().solute_moles()[1],
        t,
        2e-11,
    );
    near(report.pressure_work_into_solution_j, exact_work, 2e-7);
    near(energy(&pair, b) - before, expected_energy, 2e-7);
    assert!(energy(&pair, b) < before);
    let mut shifted = nonideal(b, 1e24);
    let shifted_report = shifted.step(0.5, step_controls(), || false).unwrap();
    assert_eq!(snapshot(&shifted), snapshot(&pair));
    assert_eq!(
        shifted_report.solvent_transfer_m3.to_bits(),
        report.solvent_transfer_m3.to_bits()
    );
    assert_eq!(
        shifted_report.pressure_work_into_solution_j.to_bits(),
        report.pressure_work_into_solution_j.to_bits()
    );
}

#[test]
fn adaptive_nonideal_window_rolls_back_late_cancellation_and_refuses_foreign_model() {
    let seed = nonideal(0.8, 0.0);
    let mut successful = seed.clone();
    let mut workspace = seed.adaptive_workspace(1 << 20).unwrap();
    let controls = DonnanAdvanceControls {
        maximum_step_s: 0.1,
        relative_tolerance: 2e-3,
        absolute_amount_tolerance_mol: 1e-11,
        absolute_water_tolerance_m3: 1e-12,
        step: step_controls(),
        ..DonnanAdvanceControls::default()
    };
    let mut polls = 0usize;
    let expected = successful
        .advance(0.4, controls, &mut workspace, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(expected.accepted_steps > 1);
    let mut pair = seed;
    let before = snapshot(&pair);
    let mut calls = 0usize;
    let error = pair
        .advance(0.4, controls, &mut workspace, || {
            calls += 1;
            calls == polls
        })
        .unwrap_err();
    assert!(matches!(
        error,
        DonnanWaterError::Network(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(calls, polls);
    assert_eq!(snapshot(&pair), before);
    let actual = pair
        .advance(0.4, controls, &mut workspace, || false)
        .unwrap();
    assert_eq!(snapshot(&pair), snapshot(&successful));
    assert_eq!(
        actual.pressure_work_into_solution_j.to_bits(),
        expected.pressure_work_into_solution_j.to_bits()
    );
    near(actual.water_balance_residual_m3, 0.0, 2e-14);
    for residual in &actual.solute_balance_residual_mol[..2] {
        near(*residual, 0.0, 2e-14);
    }
    let other = nonideal(0.4, 0.0);
    let mut foreign = other.adaptive_workspace(1 << 20).unwrap();
    let before = snapshot(&pair);
    assert!(matches!(
        pair.advance(0.1, controls, &mut foreign, || false),
        Err(DonnanWaterError::Network(OsmoticNetworkError::Invalid(_)))
    ));
    assert_eq!(snapshot(&pair), before);
}

#[test]
fn preparation_rejects_ideal_substitution_and_mismatched_original_virial_laws() {
    let source = model(0.8);
    let mut mismatched = cells(Some(&source), 0.0);
    mismatched[1] = cells(Some(&model(0.4)), 0.0)[1].clone();
    for owners in [mismatched, cells(None, 0.0)] {
        let result = DonnanOsmoticPair::new_nonideal(
            owners,
            VirialDonnan::new(vec![1, -1], Arc::clone(&source), 1e-11).unwrap(),
            [-1.0, 0.0],
            [0.1, 0.0],
            1.0,
            HydraulicPermeance::new(1.0 / RT).unwrap(),
            1 << 20,
        );
        assert!(matches!(
            result,
            Err(DonnanWaterError::Network(OsmoticNetworkError::Invalid(_)))
        ));
    }
    assert!(
        DonnanOsmoticPair::new(
            cells(Some(&source), 0.0),
            IdealDonnan::new(vec![1, -1]).unwrap(),
            [-1.0, 0.0],
            [0.1, 0.0],
            1.0,
            HydraulicPermeance::new(1.0 / RT).unwrap(),
            1 << 20
        )
        .is_err()
    );
}
