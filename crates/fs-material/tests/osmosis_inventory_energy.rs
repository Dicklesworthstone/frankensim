//! Independent finite-inventory checks of the original isothermal potential.
use fs_material::osmosis::nonideal::{NonidealSolutionError, VirialSolutionModel};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const T: f64 = 298.15;
const B: [f64; 4] = [0.002, -0.0003, -0.0003, 0.001];

fn model(reference: f64) -> VirialSolutionModel {
    VirialSolutionModel::new(T, 18e-6, reference, vec![200.0; 2], B.to_vec()).unwrap()
}

// Deliberately direct expression for moderate, separated states. No production
// difference routine or chemical-potential helper supplies this reference.
fn potential(volume: f64, amounts: &[f64; 2]) -> f64 {
    let ideal: f64 = amounts
        .iter()
        .map(|&n| {
            if n == 0.0 {
                0.0
            } else {
                n * ((n / volume).ln() - 1.0)
            }
        })
        .sum();
    let mut excess = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            excess += 0.5 * amounts[i] * B[2 * i + j] * amounts[j] / volume;
        }
    }
    R * T * (ideal + excess)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_variable_volume_and_amounts_follow_one_reference_independent_potential() {
    let (va, vb) = (0.002, 0.003);
    let (a, b) = ([0.04, 0.10], [0.12, 0.06]);
    let expected = potential(vb, &b) - potential(va, &a);
    let mut results = Vec::new();
    for reference in [1e-250, 1.0, 1e250] {
        let law = model(reference);
        let change = law
            .reference_cancelled_inventory_free_energy_change_j(va, &a, vb, &b)
            .unwrap();
        close(change, expected, 3e-12 * expected.abs());
        let reverse = law
            .reference_cancelled_inventory_free_energy_change_j(vb, &b, va, &a)
            .unwrap();
        close(change + reverse, 0.0, 2e-12);
        results.push(change);
    }
    assert_eq!(results[0], results[1]);
    assert_eq!(results[1], results[2]);

    let law = model(1.0);
    let c_before = a.map(|n| n / va);
    let c_after = b.map(|n| n / va);
    close(
        law.reference_cancelled_inventory_free_energy_change_j(va, &a, va, &b)
            .unwrap(),
        va * law
            .reference_cancelled_free_energy_change_j_m3(&c_before, &c_after)
            .unwrap(),
        2e-12,
    );
}

#[test]
fn g1_dilution_retains_osmotic_work_when_rounded_concentrations_do_not_change() {
    let volume = 0.00175_f64;
    let next_volume = f64::from_bits(volume.to_bits() + 1);
    let amount = 0.07;
    assert_eq!(amount / volume, amount / next_volume);
    let coefficient = 0.003;
    let law = VirialSolutionModel::new(T, 18e-6, 1e200, vec![200.0], vec![coefficient]).unwrap();
    let dv = next_volume - volume;
    let expected = R
        * T
        * (-amount * (dv / volume).ln_1p()
            - 0.5 * coefficient * amount * (amount / volume) * (dv / next_volume));
    let actual = law
        .reference_cancelled_inventory_free_energy_change_j(
            volume,
            &[amount],
            next_volume,
            &[amount],
        )
        .unwrap();
    assert!(actual < 0.0);
    close(actual, expected, 4e-15 * expected.abs());
    // The constitutive pressure is the exact fixed-amount volume derivative.
    let pressure = law
        .evaluate_inventories(volume, &[amount])
        .unwrap()
        .osmotic_pressure_pa;
    close(actual / dv, -pressure, 1e-9 * pressure);
}

#[test]
fn g1_inventory_derivative_agrees_with_original_chemical_potentials() {
    let law = model(1.0);
    let volume = 0.001;
    let amounts = [0.06, 0.09];
    let mut chemical = [0.0; 2];
    law.chemical_potentials_into(&amounts.map(|n| n / volume), &mut chemical)
        .unwrap();
    let h = 1e-7;
    for i in 0..2 {
        let mut minus = amounts;
        let mut plus = amounts;
        minus[i] -= h;
        plus[i] += h;
        let change = law
            .reference_cancelled_inventory_free_energy_change_j(volume, &minus, volume, &plus)
            .unwrap();
        close(change / (plus[i] - minus[i]), chemical[i], 2e-8);
    }
}

#[test]
fn g4_exact_absent_species_limits_and_original_inventory_domains_are_preserved() {
    let law = model(1.0);
    let (va, vb) = (0.001, 0.003);
    let (a, b) = ([0.0, 0.1], [0.15, 0.0]);
    close(
        law.reference_cancelled_inventory_free_energy_change_j(va, &a, vb, &b)
            .unwrap(),
        potential(vb, &b) - potential(va, &a),
        2e-12,
    );
    assert_eq!(
        law.reference_cancelled_inventory_free_energy_change_j(va, &[0.0; 2], vb, &[0.0; 2])
            .unwrap(),
        0.0
    );
    assert!(matches!(
        law.reference_cancelled_inventory_free_energy_change_j(va, &a, vb, &[1.0, 0.0]),
        Err(NonidealSolutionError::OutsideDomain { .. })
    ));
    assert!(matches!(
        law.reference_cancelled_inventory_free_energy_change_j(0.0, &a, vb, &b),
        Err(NonidealSolutionError::Invalid(_))
    ));
    assert!(matches!(
        law.reference_cancelled_inventory_free_energy_change_j(va, &a, vb, &[0.1]),
        Err(NonidealSolutionError::Shape { .. })
    ));
}
