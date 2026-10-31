//! The original reference receipts and the physical temperature snapshot have
//! distinct jobs: only the latter supplies native q knots, energy and tangent.
use super::*;
use fs_couple::material_circuit::{MaterialCircuit, motion::dielectric::ReversibleDielectricGap};

#[test]
fn native_curve_uses_actual_temperature_scaled_displacement_at_every_source_knot() {
    let card = fixture_scaled(1.0, false, 0.0, [250.0, 500.0], 0.001, 0.0);
    let (_, _, law) = initial_scaled(&card, 2.5, None, false);
    let resolved = law.at_temperature(320.0).unwrap();
    let capacitor =
        MaterialDielectricCapacitor::new_thermal(0, resolved.clone(), 0.2, 0.05).unwrap();
    let scale = 1.02;
    let reference_knots = [(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)];
    for ((field, answer), (expected_field, reference_d)) in capacitor
        .material()
        .source_path()
        .knots
        .iter()
        .zip(reference_knots)
    {
        assert_eq!(*field, expected_field);
        assert_eq!(answer.evidence.value.value, reference_d);
        for sign in [-1.0, 1.0] {
            let voltage = sign * 0.05 * field;
            let charge = sign * 0.2 * scale * reference_d;
            near(capacitor.charge_c(voltage).unwrap(), charge, 5e-16);
            near(capacitor.voltage_v(charge).unwrap(), voltage, 5e-16);
        }
    }
    near(capacitor.maximum_charge_c(), 0.2 * scale * 7.0, 5e-16);
    near(capacitor.maximum_voltage_v(), 0.2, 1e-16);
    for sign in [-1.0, 1.0] {
        let charge = sign * 0.2 * scale * 2.5;
        near(capacitor.voltage_v(charge).unwrap(), sign * 0.1, 2e-16);
        // Reference middle branch: E=2, w=1.75, epsilon_diff=0.5.
        near(
            capacitor.stored_energy_j(charge).unwrap(),
            0.01 * scale * 1.75,
            2e-16,
        );
        near(
            capacitor
                .curve()
                .differential_capacitance_f(charge)
                .unwrap(),
            4.0 * scale * 0.5,
            2e-14,
        );
        let physical = resolved.state(charge.abs() / 0.2, 0.01).unwrap();
        near(
            capacitor.stored_energy_j(charge).unwrap(),
            physical.free_energy_j,
            2e-16,
        );
        assert_ne!(physical.internal_energy_j, physical.free_energy_j);
    }
}

#[test]
fn temperature_refresh_retains_authority_and_closed_owners_cannot_drop_entropy() {
    let card = fixture_scaled(1.0, false, 0.0, [250.0, 500.0], 0.001, 0.0);
    let (_, _, law) = initial_scaled(&card, 2.5, None, false);
    let cold = law.at_temperature(300.0).unwrap();
    let hot = law.at_temperature(320.0).unwrap();
    let before = MaterialDielectricCapacitor::new_thermal(0, cold.clone(), 1.0, 1.0).unwrap();
    let after = before.at_resolved_temperature(&hot).unwrap();
    near(before.charge_c(2.0).unwrap(), 2.5, 1e-15);
    near(after.charge_c(2.0).unwrap(), 2.55, 1e-15);
    assert_eq!(after.thermal_material(), Some(&hot));
    assert_eq!(before.thermal_material(), Some(&cold));
    assert!(!std::ptr::eq(before.curve().knots(), after.curve().knots()));
    let same = after.at_resolved_temperature(&hot).unwrap();
    assert!(std::ptr::eq(same.curve().knots(), after.curve().knots()));

    // Even the scale-one reference snapshot has nonzero thermal entropy.
    assert!(MaterialDielectricCapacitor::new(0, cold.dielectric().clone(), 1.0, 1.0).is_err());
    assert!(ReversibleDielectricGap::new(cold.dielectric().clone(), 0, 1.0, [0.5, 1.5]).is_err());
    let field_point = point(320.0)
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    assert!(before.at_temperature_query_point(&field_point).is_err());
    assert!(
        MaterialCircuit::new_with_dielectric_capacitors(
            CircuitGraph {
                node_count: 2,
                branches: vec![(1, 0, Branch::Capacitor { farads: f64::NAN })],
                transformers: vec![],
            },
            vec![],
            vec![],
            vec![],
            vec![],
            vec![before],
            vec![],
            vec![2.5],
            &[],
            || false,
        )
        .is_err()
    );
}
