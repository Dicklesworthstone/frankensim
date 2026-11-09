//! Analytical checks of the exact production conservation core. The two-cell
//! elimination below is a test oracle, not a replacement for the native LU.
use super::*;

fn checkpoint() -> core::result::Result<(), ()> { Ok(()) }
fn controls() -> IsolatedThermalImpedanceControls {
    IsolatedThermalImpedanceControls {
        cell_power_tolerance_w: 1e-12,
        component_power_tolerance_w: 1e-12,
        component_energy_tolerance_j: 1e-12,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}
fn pair(slopes: [f64; 2], conductance: f64) -> (ConservedThermalEnergy, [f64; 4]) {
    let tangent = [
        conductance * slopes[0], -conductance * slopes[1],
        -conductance * slopes[0], conductance * slopes[1],
    ];
    let energy = ConservedThermalEnergy::new(
        vec![vec![0, 1]], &slopes, &tangent, controls(), &mut checkpoint,
    ).unwrap();
    (energy, tangent)
}
fn solve_pair(matrix: &ConservedThermalMatrix, drive: &ConservedThermalDrive) -> [C64; 2] {
    assert_eq!(matrix.values[0], C64::ONE);
    assert_eq!(matrix.values[1], C64::ONE);
    assert_eq!(matrix.row_scales[0], 1.0);
    let total = drive.rhs[0];
    let c = matrix.values[2];
    let d = matrix.values[3];
    let rhs = drive.rhs[1].scale(1.0 / matrix.row_scales[1]);
    let second = (rhs - c * total) / (d - c);
    [total - second, second]
}
fn original_rows(tangent: &[f64; 4], omega: f64, heat: [C64; 2], energy: [C64; 2]) {
    for row in 0..2 {
        let residual = energy[0].scale(tangent[2 * row])
            + energy[1].scale(tangent[2 * row + 1])
            + C64::new(omega * energy[row].im, -omega * energy[row].re)
            - heat[row];
        near(residual, C64::ZERO, 5e-14);
    }
}

#[test]
fn dc_retains_extensive_energy_with_unequal_heat_capacities() {
    // C0=2, C1=3 J/K; G=4 W/K. Fixing one temperature would be wrong.
    let (owner, tangent) = pair([0.5, 1.0 / 3.0], 4.0);
    let matrix = owner.matrix(&tangent, 0.0, &mut checkpoint).unwrap();
    for total in [0.0, 10.0] {
        let heat = [C64::ONE, -C64::ONE];
        let inventory = [C64::from_re(total)];
        let drive = owner.drive(0.0, &heat, Some(&inventory), &mut checkpoint).unwrap();
        let energy = solve_pair(&matrix, &drive);
        near(energy[0], C64::from_re(0.4 * total + 0.3), 2e-15);
        near(energy[1], C64::from_re(0.6 * total - 0.3), 2e-15);
        near(energy[0].scale(0.5) - energy[1].scale(1.0 / 3.0), C64::from_re(0.25), 1e-15);
        let report = owner.admit(&drive, &energy, &mut checkpoint).unwrap();
        near(report[0].energy_j, C64::from_re(total), 1e-15);
        original_rows(&tangent, 0.0, heat, energy);
    }
    let drive = owner.drive(0.0, &[C64::ZERO; 2], Some(&[C64::from_re(10.0)]), &mut checkpoint).unwrap();
    let energy = solve_pair(&matrix, &drive);
    near(energy[0], C64::from_re(4.0), 1e-15);
    near(energy[1], C64::from_re(6.0), 1e-15);
    owner.admit(&drive, &energy, &mut checkpoint).unwrap();
}

#[test]
fn balanced_harmonic_heat_has_the_same_conserved_limit_at_tiny_frequencies() {
    let (owner, tangent) = pair([0.5, 1.0 / 3.0], 4.0);
    let heat = [C64::ONE, -C64::ONE];
    for omega in [1.0, 1e-12, 1e-100] {
        let matrix = owner.matrix(&tangent, omega, &mut checkpoint).unwrap();
        let drive = owner.drive(omega, &heat, None, &mut checkpoint).unwrap();
        let energy = solve_pair(&matrix, &drive);
        let expected = C64::new(10.0 / 3.0, -omega).recip();
        near(energy[0], expected, 2e-16);
        near(energy[1], -expected, 2e-16);
        let report = owner.admit(&drive, &energy, &mut checkpoint).unwrap();
        assert_eq!(report[0].energy_j, C64::ZERO);
        assert_eq!(report[0].energy_residual_j, C64::ZERO);
        original_rows(&tangent, omega, heat, energy);
    }
}

#[test]
fn net_harmonic_heating_retains_the_physical_storage_pole() {
    let (owner, tangent) = pair([0.5, 1.0 / 3.0], 4.0);
    let heat = [C64::ONE, C64::ZERO];
    for omega in [2.0, 0.25] {
        let matrix = owner.matrix(&tangent, omega, &mut checkpoint).unwrap();
        let drive = owner.drive(omega, &heat, None, &mut checkpoint).unwrap();
        let energy = solve_pair(&matrix, &drive);
        let total = C64::new(0.0, 1.0 / omega);
        let first = C64::new(4.0 / 3.0, -omega) / C64::new(10.0 / 3.0, -omega) * total;
        near(energy[0], first, 1e-15);
        near(energy[1], total - first, 1e-15);
        let report = owner.admit(&drive, &energy, &mut checkpoint).unwrap();
        near(report[0].energy_j, total, 1e-15);
        near(report[0].net_heat_w, C64::ONE, 0.0);
        original_rows(&tangent, omega, heat, energy);
        assert!(owner.drive(omega, &heat, Some(&[C64::ZERO]), &mut checkpoint).is_err());
    }
    assert!(matches!(owner.drive(0.0, &heat, None, &mut checkpoint),
        Err(ThermalEnergyConstraintError::IncompatibleHeat { component: 0, .. })));
}

#[test]
fn singleton_latent_and_sensible_cells_keep_storage_without_a_bath() {
    for slope in [0.0, 0.2] {
        let owner = ConservedThermalEnergy::new(
            vec![vec![0]], &[slope], &[0.0], controls(), &mut checkpoint,
        ).unwrap();
        for omega in [0.0, 2.0] {
            let matrix = owner.matrix(&[0.0], omega, &mut checkpoint).unwrap();
            assert_eq!(matrix.values, [C64::ONE]);
            let heat = if omega == 0.0 { C64::ZERO } else { C64::new(3.0, 4.0) };
            let prescribed = [C64::from_re(7.0)];
            let drive = owner.drive(omega, &[heat], (omega == 0.0).then_some(prescribed.as_slice()), &mut checkpoint).unwrap();
            let expected = if omega == 0.0 { prescribed[0] } else { C64::new(-2.0, 1.5) };
            near(drive.rhs[0], expected, 0.0);
            owner.admit(&drive, &[expected], &mut checkpoint).unwrap();
        }
    }
}

#[test]
fn one_latent_parcel_is_fixed_by_inventory_but_several_are_not() {
    let (owner, tangent) = pair([0.0, 0.5], 4.0);
    let heat = [C64::ONE, -C64::ONE];
    let matrix = owner.matrix(&tangent, 0.0, &mut checkpoint).unwrap();
    let drive = owner.drive(0.0, &heat, None, &mut checkpoint).unwrap();
    let energy = solve_pair(&matrix, &drive);
    near(energy[0], C64::from_re(0.5), 0.0);
    near(energy[1], C64::from_re(-0.5), 0.0);
    owner.admit(&drive, &energy, &mut checkpoint).unwrap();
    original_rows(&tangent, 0.0, heat, energy);

    let (owner, tangent) = pair([0.0, 0.0], 4.0);
    assert!(matches!(owner.matrix(&tangent, 0.0, &mut checkpoint),
        Err(ThermalEnergyConstraintError::Underdetermined { component: 0, latent_cells: 2 })));
    let matrix = owner.matrix(&tangent, 2.0, &mut checkpoint).unwrap();
    let drive = owner.drive(2.0, &heat, None, &mut checkpoint).unwrap();
    let energy = solve_pair(&matrix, &drive);
    near(energy[0], C64::new(0.0, 0.5), 0.0);
    near(energy[1], C64::new(0.0, -0.5), 0.0);
    owner.admit(&drive, &energy, &mut checkpoint).unwrap();
    original_rows(&tangent, 2.0, heat, energy);
}

#[test]
fn disconnected_components_cannot_cancel_each_others_net_heat() {
    let owner = ConservedThermalEnergy::new(
        vec![vec![0], vec![1]], &[1.0; 2], &[0.0; 4], controls(), &mut checkpoint,
    ).unwrap();
    assert!(matches!(owner.drive(0.0, &[C64::ONE, -C64::ONE], None, &mut checkpoint),
        Err(ThermalEnergyConstraintError::IncompatibleHeat { component: 0, .. })));
    let inventories = [C64::from_re(3.0), C64::from_re(-4.0)];
    let drive = owner.drive(0.0, &[C64::ZERO; 2], Some(&inventories), &mut checkpoint).unwrap();
    assert_eq!(drive.rhs, inventories);
    let report = owner.admit(&drive, &inventories, &mut checkpoint).unwrap();
    assert_eq!(report.len(), 2);
    assert_eq!(report[0].energy_j, inventories[0]);
    assert_eq!(report[1].energy_j, inventories[1]);
}

#[test]
fn compensated_original_heat_sum_and_unrounded_energy_target_are_authoritative() {
    let owner = ConservedThermalEnergy::new(
        vec![vec![0, 1, 2]], &[1.0; 3],
        &[1.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 1.0],
        controls(), &mut checkpoint,
    ).unwrap();
    let heat = [C64::from_re(1e16), C64::ONE, C64::from_re(-1e16)];
    assert!(matches!(owner.drive(0.0, &heat, None, &mut checkpoint),
        Err(ThermalEnergyConstraintError::IncompatibleHeat { net_power_w: C64 { re: 1.0, im: 0.0 }, .. })));
    let drive = owner.drive(2.0, &heat, None, &mut checkpoint).unwrap();
    assert_eq!(drive.rhs[0], C64::new(0.0, 0.5));

    let (owner, _) = pair([1.0; 2], 1.0);
    let drive = owner.drive(1.0, &[C64::from_re(1e16), C64::ONE], None, &mut checkpoint).unwrap();
    assert_eq!(drive.rhs[0], C64::new(0.0, 1e16));
    assert!(matches!(owner.admit(&drive, &[C64::new(0.0, 1e16), C64::ZERO], &mut checkpoint),
        Err(ThermalEnergyConstraintError::EnergyResidual { residual_j: C64 { re: 0.0, im: -1.0 }, .. })));
}

#[test]
fn absolute_energy_and_power_gates_are_independent() {
    let (owner, _) = pair([1.0; 2], 1.0);
    let dc = owner.drive(0.0, &[C64::ZERO; 2], None, &mut checkpoint).unwrap();
    assert!(matches!(owner.admit(&dc, &[C64::ONE; 2], &mut checkpoint),
        Err(ThermalEnergyConstraintError::EnergyResidual { .. })));
    let ac = owner.drive(1e8, &[C64::ZERO; 2], None, &mut checkpoint).unwrap();
    assert!(matches!(owner.admit(&ac, &[C64::from_re(1e-13), C64::ZERO], &mut checkpoint),
        Err(ThermalEnergyConstraintError::PowerResidual { .. })));
}

#[test]
fn original_partition_and_conservative_tangent_cannot_be_repaired_silently() {
    assert!(ConservedThermalEnergy::new(
        vec![vec![0, 0], vec![1]], &[1.0; 2], &[0.0; 4], controls(), &mut checkpoint,
    ).is_err());
    assert!(ConservedThermalEnergy::new(
        vec![vec![0]], &[1.0; 2], &[0.0; 4], controls(), &mut checkpoint,
    ).is_err());
    assert!(ConservedThermalEnergy::new(
        vec![vec![0], vec![1]], &[1.0; 2], &[1.0, -1.0, -1.0, 1.0], controls(), &mut checkpoint,
    ).is_err());
    assert!(matches!(ConservedThermalEnergy::new(
        vec![vec![0, 1]], &[1.0; 2], &[1.0, 0.0, 0.0, 0.0], controls(), &mut checkpoint,
    ), Err(ThermalEnergyConstraintError::Nonconservative { .. })));
}

#[test]
fn nonfinite_overflow_underflow_and_cancellation_refuse_without_a_response() {
    let (owner, tangent) = pair([1.0; 2], 1.0);
    assert!(owner.matrix(&tangent, -1.0, &mut checkpoint).is_err());
    assert!(owner.drive(1.0, &[C64::new(f64::NAN, 0.0); 2], None, &mut checkpoint).is_err());
    assert!(owner.drive(1e-300, &[C64::from_re(1e300), C64::ZERO], None, &mut checkpoint).is_err());
    assert!(owner.drive(2.0, &[C64::from_re(f64::from_bits(1)), C64::ZERO], None, &mut checkpoint).is_err());
    assert!(matches!(owner.matrix(&tangent, 1.0, &mut || Err("stop")),
        Err(ThermalEnergyConstraintError::Checkpoint("stop"))));
    assert!(matches!(owner.drive(1.0, &[C64::ZERO; 2], None, &mut || Err("stop")),
        Err(ThermalEnergyConstraintError::Checkpoint("stop"))));
    let drive = owner.drive(1.0, &[C64::ZERO; 2], None, &mut checkpoint).unwrap();
    assert!(matches!(owner.admit(&drive, &[C64::ZERO; 2], &mut || Err("stop")),
        Err(ThermalEnergyConstraintError::Checkpoint("stop"))));
    assert!(owner.admit(&drive, &[C64::ZERO; 2], &mut checkpoint).is_ok());
}
