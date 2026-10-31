use super::*;
use crate::osmosis::{IdealDiluteSolution, OsmoticPressureVolume};
use fs_material::osmosis::donnan::DonnanPotential;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const T: f64 = 298.15;
const RT: f64 = R * T;

fn cell(volume: f64, amounts: &[f64], pressure: f64, elastance: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(T, 18e-6, volume, amounts.to_vec()).unwrap(),
        OsmoticPressureVolume::new(volume, pressure, elastance).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}
fn pair(
    cells: [OsmoticCompartment; 2],
    charges: Vec<i32>,
    fixed: [f64; 2],
    trapped: [f64; 2],
    lp: f64,
) -> DonnanOsmoticPair {
    DonnanOsmoticPair::new(
        cells,
        IdealDonnan::new(charges).unwrap(),
        fixed,
        trapped,
        1.0,
        HydraulicPermeance::new(lp).unwrap(),
        1 << 20,
    )
    .unwrap()
}
fn reference() -> DonnanOsmoticPair {
    pair(
        [
            cell(1.0, &[1.0, 0.0], 0.0, RT),
            cell(1.0, &[1.0, 1.0], 0.0, RT),
        ],
        vec![1, -1],
        [-1.0, 0.0],
        [0.0, 0.0],
        1.0 / RT,
    )
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn snapshot(system: &DonnanOsmoticPair) -> Vec<u64> {
    let mut result = vec![system.time_s().to_bits()];
    for cell in system.compartments() {
        result.push(cell.solution().solvent_volume_m3().to_bits());
        result.extend(cell.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    let equilibrium = system.equilibrium();
    result.extend(equilibrium.osmotic_pressure_pa.map(f64::to_bits));
    result.extend(equilibrium.charge_residual_mol.map(f64::to_bits));
    match equilibrium.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            voltage_left_minus_right_v,
        } => result.extend([
            0,
            dimensionless_left_minus_right.to_bits(),
            voltage_left_minus_right_v.to_bits(),
        ]),
        DonnanPotential::NegativeInfinite => result.push(1),
        DonnanPotential::PositiveInfinite => result.push(2),
        DonnanPotential::Indeterminate => result.push(3),
    }
    result
}
fn energy(system: &DonnanOsmoticPair) -> f64 {
    system
        .compartments()
        .iter()
        .enumerate()
        .map(|(side, cell)| {
            let volume = cell.solution().solvent_volume_m3();
            let delta = volume - cell.pressure.reference_volume_m3;
            let mechanical = cell.pressure.reference_pressure_pa * delta
                + 0.5 * cell.pressure.elastance_pa_per_m3 * delta * delta;
            let ideal = |amount: f64| {
                if amount == 0.0 {
                    0.0
                } else {
                    amount * ((amount / volume).ln() - 1.0)
                }
            };
            mechanical / RT
                + cell
                    .solution()
                    .solute_moles()
                    .iter()
                    .map(|&n| ideal(n))
                    .sum::<f64>()
                + ideal(system.trapped_particles_mol()[side])
        })
        .sum()
}

#[test]
fn preparation_equilibrates_finite_ions_without_changing_water_or_species_totals() {
    let system = reference();
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    assert_eq!(
        (left.solvent_volume_m3(), right.solvent_volume_m3()),
        (1.0, 1.0)
    );
    assert_eq!(system.mobile_totals_mol(), &[2.0, 1.0]);
    near(left.solute_moles()[0], 4.0 / 3.0, 3e-13);
    near(left.solute_moles()[1], 1.0 / 3.0, 3e-13);
    near(right.solute_moles()[0], 2.0 / 3.0, 3e-13);
    near(right.solute_moles()[1], 2.0 / 3.0, 3e-13);
    assert_eq!(system.time_s(), 0.0);
}

#[test]
fn pure_water_linear_compliance_has_exact_backward_euler_step_and_mechanical_work() {
    let cells = [cell(1.0, &[], 2.0, 3.0), cell(1.0, &[], 0.0, 2.0)];
    let mut system = DonnanOsmoticPair::new(
        cells,
        IdealDonnan::new(vec![]).unwrap(),
        [0.0; 2],
        [0.0; 2],
        2.0,
        HydraulicPermeance::new(1.0).unwrap(),
        1 << 20,
    )
    .unwrap();
    let controls = DonnanStepControls {
        maximum_root_evaluations: 2,
        ..DonnanStepControls::default()
    };
    let report = system.step(0.1, controls, || false).unwrap();
    // x=dt*A*Lp*dp/(1+dt*A*Lp*(KL+KR))=0.2.
    near(report.solvent_transfer_m3, 0.2, 2e-14);
    near(report.pressure_work_into_solution_j, 0.3, 2e-14);
    near(
        system.compartments()[0].solution().solvent_volume_m3(),
        0.8,
        2e-14,
    );
    assert_eq!(report.root_evaluations, 2);
    assert_eq!(report.ion_root_evaluations, 0);
    assert_eq!(report.equilibrium.potential, DonnanPotential::Indeterminate);
}

#[test]
fn coupled_ionic_water_equilibrium_matches_a_known_finite_bath_partition() {
    let mut system = pair(
        [
            cell(0.8, &[1.0, 0.0], RT / 3.0, 0.0),
            cell(1.2, &[1.0, 1.0], 0.0, 0.0),
        ],
        vec![1, -1],
        [-1.0, 0.0],
        [0.0; 2],
        1.0 / RT,
    );
    for _ in 0..80 {
        system
            .step(1.0, DonnanStepControls::default(), || false)
            .unwrap();
    }
    near(
        system.compartments()[0].solution().solvent_volume_m3(),
        1.0,
        2e-8,
    );
    near(
        system.compartments()[1].solution().solvent_volume_m3(),
        1.0,
        2e-8,
    );
    near(
        system.equilibrium().osmotic_pressure_pa[0] - system.equilibrium().osmotic_pressure_pa[1],
        RT / 3.0,
        3e-5,
    );
    assert_eq!(system.fixed_charge_mol(), [-1.0, 0.0]);
}

#[test]
fn very_large_hydraulic_permeance_uses_water_bracket_accuracy_without_overshoot() {
    let mut system = pair(
        [
            cell(0.8, &[1.0, 0.0], RT / 3.0, 0.0),
            cell(1.2, &[1.0, 1.0], 0.0, 0.0),
        ],
        vec![1, -1],
        [-1.0, 0.0],
        [0.0; 2],
        1e10 / RT,
    );
    let report = system
        .step(1.0, DonnanStepControls::default(), || false)
        .unwrap();
    near(
        system.compartments()[0].solution().solvent_volume_m3(),
        1.0,
        2e-8,
    );
    assert!(report.solvent_transfer_m3 < 0.0);
    assert!(report.root_evaluations <= 128);
    near(report.water_balance_residual_m3, 0.0, 1e-14);
    assert!(
        report
            .equilibrium
            .charge_residual_mol
            .iter()
            .all(|q| q.abs() < 2e-12)
    );
}

#[test]
fn counterion_only_swelling_has_an_exact_water_root_and_distinct_trapped_pressure() {
    let mut system = pair(
        [cell(1.0, &[1.0], 0.0, RT), cell(1.0, &[0.0], 0.0, RT)],
        vec![1],
        [-1.0, 0.0],
        [0.1, 0.0],
        0.1 / RT,
    );
    let before = energy(&system);
    let report = system
        .step(0.5, DonnanStepControls::default(), || false)
        .unwrap();
    let expected = (1.0 - 1.2_f64.sqrt()) / 2.0;
    near(report.solvent_transfer_m3, expected, 2e-12);
    assert_eq!(
        report.equilibrium.potential,
        DonnanPotential::NegativeInfinite
    );
    assert_eq!(system.compartments()[0].solution().solute_moles(), &[1.0]);
    assert_eq!(system.compartments()[1].solution().solute_moles(), &[0.0]);
    let left = system.compartments()[0].solution();
    near(
        report.equilibrium.osmotic_pressure_pa[0] - left.osmotic_pressure_pa().unwrap(),
        RT * 0.1 / left.solvent_volume_m3(),
        2e-9,
    );
    assert!(report.pressure_work_into_solution_j < 0.0);
    assert!(energy(&system) < before);
}

#[test]
fn multivalent_steps_conserve_every_inventory_charge_and_decrease_free_energy() {
    let cells = [
        cell(1.0, &[1.0, 1.0, 0.5], 0.2 * RT, 0.4 * RT),
        cell(1.0, &[2.0, 2.0, 0.5], 0.0, 0.4 * RT),
    ];
    let mut forward = pair(
        cells.clone(),
        vec![1, -1, 2],
        [-1.0; 2],
        [0.1, 0.3],
        0.2 / RT,
    );
    let mut reverse = pair(
        [cells[1].clone(), cells[0].clone()],
        vec![1, -1, 2],
        [-1.0; 2],
        [0.3, 0.1],
        0.2 / RT,
    );
    let mut previous = energy(&forward);
    for _ in 0..20 {
        let report = forward
            .step(0.1, DonnanStepControls::default(), || false)
            .unwrap();
        reverse
            .step(0.1, DonnanStepControls::default(), || false)
            .unwrap();
        let current = energy(&forward);
        assert!(current <= previous + 1e-10 * (1.0 + previous.abs()));
        previous = current;
        near(
            forward
                .compartments()
                .iter()
                .map(|c| c.solution().solvent_volume_m3())
                .sum(),
            2.0,
            1e-13,
        );
        for (i, total) in [3.0, 3.0, 1.0].into_iter().enumerate() {
            near(
                forward
                    .compartments()
                    .iter()
                    .map(|c| c.solution().solute_moles()[i])
                    .sum(),
                total,
                2e-14,
            );
            near(report.solute_balance_residual_mol[i], 0.0, 2e-14);
        }
        for side in 0..2 {
            let a = forward.compartments()[side].solution();
            let b = reverse.compartments()[1 - side].solution();
            near(a.solvent_volume_m3(), b.solvent_volume_m3(), 3e-10);
            for (&na, &nb) in a.solute_moles().iter().zip(b.solute_moles()) {
                near(na, nb, 3e-10);
            }
            near(report.equilibrium.charge_residual_mol[side], 0.0, 1e-11);
        }
    }
    assert!((forward.compartments()[0].solution().solvent_volume_m3() - 1.0).abs() > 1e-3);
}

#[test]
fn an_unreachable_water_root_refuses_instead_of_clipping_a_volume_floor() {
    let mut cells = [cell(1.0, &[], 1.0, 0.0), cell(1.0, &[], 0.0, 0.0)];
    cells[0].minimum_volume_m3 = 0.25;
    let mut system = pair(cells, vec![], [0.0; 2], [0.0; 2], 1.0);
    let before = snapshot(&system);
    let failure = system.step(2.0, DonnanStepControls::default(), || false);
    assert!(matches!(
        failure,
        Err(DonnanWaterError::Network(OsmoticNetworkError::Unbracketed))
    ));
    assert_eq!(snapshot(&system), before);
}

#[test]
fn late_cancellation_and_both_work_limits_preserve_state_and_retry_exactly() {
    let mut system = reference();
    let before = snapshot(&system);
    let outer = DonnanStepControls {
        maximum_root_evaluations: 2,
        ..DonnanStepControls::default()
    };
    assert!(matches!(
        system.step(1.0, outer, || false),
        Err(DonnanWaterError::Network(
            OsmoticNetworkError::NonlinearConvergence
        ))
    ));
    assert_eq!(snapshot(&system), before);
    let inner = DonnanStepControls {
        ion_equilibrium: DonnanControls {
            maximum_iterations: 1,
            ..DonnanControls::default()
        },
        ..DonnanStepControls::default()
    };
    assert!(matches!(
        system.step(1.0, inner, || false),
        Err(DonnanWaterError::Ions(DonnanError::NonConvergence))
    ));
    assert_eq!(snapshot(&system), before);
    let mut uninterrupted = system.clone();
    let mut final_poll = 0;
    uninterrupted
        .step(1.0, DonnanStepControls::default(), || {
            final_poll += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        system.step(1.0, DonnanStepControls::default(), || {
            polls += 1;
            polls == final_poll
        }),
        Err(DonnanWaterError::Network(OsmoticNetworkError::Cancelled))
    ));
    assert_eq!(snapshot(&system), before);
    system
        .step(1.0, DonnanStepControls::default(), || false)
        .unwrap();
    assert_eq!(snapshot(&system), snapshot(&uninterrupted));
}

#[test]
fn preparation_counts_ion_workspace_and_rejects_incompatible_nonideal_basis() {
    use fs_material::osmosis::{AqueousSolution, nonideal::VirialSolutionModel};
    use std::sync::Arc;
    let source = reference();
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    assert!(
        source.workspace_bytes()
            >= core::mem::size_of::<DonnanOsmoticPair>()
                + law.heap_bytes()
                + 7 * 2 * core::mem::size_of::<f64>()
    );
    let refused = DonnanOsmoticPair::new(
        source.compartments().clone(),
        law.clone(),
        [-1.0, 0.0],
        [0.0; 2],
        1.0,
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        source.workspace_bytes() - 1,
    );
    assert!(matches!(
        refused,
        Err(DonnanWaterError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
    let model = Arc::new(
        VirialSolutionModel::new(T, 18e-6, 1.0, vec![10.0; 2], vec![0.1, 0.0, 0.0, 0.1]).unwrap(),
    );
    let nonideal = AqueousSolution::nonideal(model, 1.0, vec![1.0, 0.0]).unwrap();
    let left = OsmoticCompartment::new(
        nonideal,
        OsmoticPressureVolume::new(1.0, 0.0, RT).unwrap(),
        1e-6,
    )
    .unwrap();
    assert!(
        DonnanOsmoticPair::new(
            [left, cell(1.0, &[1.0, 1.0], 0.0, RT)],
            law,
            [-1.0, 0.0],
            [0.0; 2],
            1.0,
            HydraulicPermeance::new(1.0 / RT).unwrap(),
            1 << 20
        )
        .is_err()
    );
}
