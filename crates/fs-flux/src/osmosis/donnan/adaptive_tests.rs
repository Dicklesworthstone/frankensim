use super::*;
use crate::osmosis::{IdealDiluteSolution, OsmoticPressureVolume};
use fs_material::osmosis::donnan::DonnanPotential;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const T: f64 = 298.15;

fn cell(amounts: &[f64], pressure: f64, elastance: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(T, 18e-6, 1.0, amounts.to_vec()).unwrap(),
        OsmoticPressureVolume::new(1.0, pressure, elastance).unwrap(),
        1e-6,
    )
    .unwrap()
}
fn pair(amounts: &[f64], valences: Vec<i32>) -> DonnanOsmoticPair {
    DonnanOsmoticPair::new(
        [cell(amounts, 1.0, 1.0), cell(amounts, 0.0, 1.0)],
        IdealDonnan::new(valences).unwrap(),
        [0.0; 2],
        [0.0; 2],
        1.0,
        HydraulicPermeance::new(1.0).unwrap(),
        1 << 20,
    )
    .unwrap()
}
fn controls(tolerance: f64) -> DonnanAdvanceControls {
    DonnanAdvanceControls {
        relative_tolerance: tolerance,
        ..DonnanAdvanceControls::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn water(pair: &DonnanOsmoticPair) -> f64 {
    pair.compartments()[0].solution().solvent_volume_m3()
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
fn refined_compliant_relaxation_tracks_analytic_time_law_and_integrated_work() {
    // dx/dt = 1 - 2x, x(0)=0. Algebraic BE accuracy cannot guarantee this transient.
    let exact = 0.5 * (1.0 - (-2.0_f64).exp());
    let mut errors = [0.0; 2];
    let mut counts = [0; 2];
    for (i, tolerance) in [1e-2, 1e-6].into_iter().enumerate() {
        let mut system = pair(&[], vec![]);
        let mut workspace = system.adaptive_workspace(1 << 20).unwrap();
        let report = system
            .advance(1.0, controls(tolerance), &mut workspace, || false)
            .unwrap();
        errors[i] = (report.solvent_transfer_m3 - exact).abs();
        counts[i] = report.accepted_steps;
        let x = report.solvent_transfer_m3;
        near(report.pressure_work_into_solution_j, x - x * x, 2e-13);
        near(report.water_balance_residual_m3, 0.0, 2e-14);
        assert_eq!(report.interval_s, [0.0, 1.0]);
        assert_eq!(system.time_s(), 1.0);
        assert!(report.maximum_scaled_error <= 1.0);
        assert_eq!(
            report.step_evaluations,
            3 * (report.accepted_steps + report.rejected_steps)
        );
        assert!(report.root_evaluations >= report.step_evaluations);
        assert_eq!(report.ion_root_evaluations, 0);
        assert!(report.endpoint.interval_s[0] > 0.0);
        assert_eq!(report.endpoint.interval_s[1], 1.0);
    }
    assert!(errors[1] < 2e-4, "refined error={:e}", errors[1]);
    assert!(errors[1] < errors[0] / 20.0);
    assert!(counts[1] > counts[0]);
}

#[test]
fn every_mobile_amount_controls_ionic_repartition_even_with_loose_water_budget() {
    // No fixed charge: both baths retain the same concentrations as water
    // redistributes, so n_left_i=N_i*V_left/2 and the same analytic x(t) holds.
    let amounts = [1.0, 1.0, 1e-12, 1e-12, 0.0, 0.0];
    let mut system = pair(&amounts, vec![1, -1, 2, -2, 3, -3]);
    let mut workspace = system.adaptive_workspace(1 << 20).unwrap();
    let physical = DonnanAdvanceControls {
        relative_tolerance: 1e-8,
        absolute_water_tolerance_m3: 1.0,
        absolute_amount_tolerance_mol: 1e-6,
        ..DonnanAdvanceControls::default()
    };
    let report = system
        .advance(1.0, physical, &mut workspace, || false)
        .unwrap();
    let exact_volume = 0.5 * (1.0 + (-2.0_f64).exp());
    assert!(report.accepted_steps > 20);
    near(water(&system), exact_volume, 2e-4);
    for (i, &initial) in amounts.iter().enumerate() {
        let left = system.compartments()[0].solution().solute_moles()[i];
        let right = system.compartments()[1].solution().solute_moles()[i];
        near(left, initial * exact_volume, 2e-4 * initial);
        near(left + right, 2.0 * initial, 4e-15 * initial);
        near(report.solute_balance_residual_mol[i], 0.0, 4e-15 * initial);
        if initial == 0.0 {
            assert_eq!((left, right), (0.0, 0.0));
        }
    }
    assert!(report.ion_root_evaluations > 0);
}

#[test]
fn counterion_only_swelling_resolves_analytic_square_root_without_salt_floors() {
    // dV_left/dt = 0.1*(1 counterion + 0.2 trapped osmoles)/V_left.
    let mut system = DonnanOsmoticPair::new(
        [cell(&[1.0], 0.0, 0.0), cell(&[0.0], 0.0, 0.0)],
        IdealDonnan::new(vec![1]).unwrap(),
        [-1.0, 0.0],
        [0.2, 0.0],
        1.0,
        HydraulicPermeance::new(0.1 / (R * T)).unwrap(),
        1 << 20,
    )
    .unwrap();
    let mut workspace = system.adaptive_workspace(1 << 20).unwrap();
    let report = system
        .advance(1.0, controls(1e-7), &mut workspace, || false)
        .unwrap();
    near(water(&system), 1.24_f64.sqrt(), 3e-5);
    assert_eq!(system.compartments()[0].solution().solute_moles(), &[1.0]);
    assert_eq!(system.compartments()[1].solution().solute_moles(), &[0.0]);
    assert_eq!(
        system.equilibrium().potential,
        DonnanPotential::NegativeInfinite
    );
    assert_eq!(report.solute_balance_residual_mol[0], 0.0);
    assert_eq!(report.pressure_work_into_solution_j, 0.0);
}

#[test]
fn late_cancellation_and_work_exhaustion_restore_inventories_equilibrium_and_clock() {
    let charged = || {
        DonnanOsmoticPair::new(
            [cell(&[1.0, 0.0], 0.0, R * T), cell(&[1.0, 1.0], 0.0, R * T)],
            IdealDonnan::new(vec![1, -1]).unwrap(),
            [-1.0, 0.0],
            [0.0; 2],
            1.0,
            HydraulicPermeance::new(1.0 / (R * T)).unwrap(),
            1 << 20,
        )
        .unwrap()
    };
    let mut successful = charged();
    let mut workspace = successful.adaptive_workspace(1 << 20).unwrap();
    let physical = DonnanAdvanceControls {
        maximum_step_s: 0.1,
        ..controls(1e-3)
    };
    let mut polls = 0usize;
    let result = successful
        .advance(0.4, physical, &mut workspace, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(result.accepted_steps > 1);
    for fail_at_end in [true, false] {
        let mut system = charged();
        let before = snapshot(&system);
        let mut calls = 0;
        let config = if fail_at_end {
            physical
        } else {
            DonnanAdvanceControls {
                maximum_steps: 1,
                ..physical
            }
        };
        let error = system
            .advance(0.4, config, &mut workspace, || {
                calls += 1;
                fail_at_end && calls == polls
            })
            .unwrap_err();
        assert!(matches!(
            error,
            DonnanWaterError::Network(
                OsmoticNetworkError::Cancelled | OsmoticNetworkError::WorkBudget
            )
        ));
        assert_eq!(snapshot(&system), before);
        // A refused window does not poison the prepared workspace.
        system
            .advance(0.4, physical, &mut workspace, || false)
            .unwrap();
        assert_eq!(snapshot(&system), snapshot(&successful));
    }
}

#[test]
fn temporal_floor_workspace_budget_and_shape_refuse_without_publication() {
    let mut system = pair(&[], vec![]);
    let before = snapshot(&system);
    let mut workspace = system.adaptive_workspace(1 << 20).unwrap();
    assert!(matches!(
        system.adaptive_workspace(workspace.workspace_bytes() - 1),
        Err(DonnanWaterError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
    let physical = DonnanAdvanceControls {
        minimum_step_s: 1.0,
        ..controls(1e-10)
    };
    assert!(matches!(
        system.advance(1.0, physical, &mut workspace, || false),
        Err(DonnanWaterError::TemporalAccuracy { .. })
    ));
    assert_eq!(snapshot(&system), before);
    let ionic = pair(&[1.0, 1.0], vec![1, -1]);
    let mut wrong_shape = ionic.adaptive_workspace(1 << 20).unwrap();
    assert!(
        system
            .advance(1.0, controls(1e-3), &mut wrong_shape, || false)
            .is_err()
    );
    assert_eq!(snapshot(&system), before);
}
