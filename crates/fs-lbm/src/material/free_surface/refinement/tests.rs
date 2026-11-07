use super::*;

fn report() -> DyadicRefinementReport3 {
    let before = DyadicRefinementInventory3 {
        mass: 3.0,
        first_mass_moment: [2.0, 3.0, 4.0],
        momentum: [0.2, -0.3, 0.4],
        angular_momentum: [-0.4, 0.5, 0.6],
        kinetic_energy: 0.07,
    };
    DyadicRefinementReport3 {
        source_cells: [4; 3],
        target_cells: [8; 3],
        source_active_tiles: 1,
        target_active_tiles: 8,
        before,
        after: before,
        absolute_mass_defect: 2e-15,
        absolute_first_moment_defect: 3e-15,
        absolute_momentum_defect: 4e-15,
        absolute_angular_momentum_defect: 5e-15,
        kinetic_energy_change: -6e-15,
        absolute_kinetic_energy_defect: 6e-15,
        maximum_viscous_stress_defect: 7e-15,
        source_geometry_epoch: 9,
        target_geometry_epoch: 10,
    }
}

#[test]
fn physical_receipt_uses_original_mass_length_velocity_and_pressure_units() {
    let native = report();
    let units = [2.5, 7.5, 10.0, 30.0, 40.0, 16.0];
    let physical = receipt(native, units, [1e-10; 6]).unwrap();
    assert_eq!(physical.before.mass_kg, 7.5);
    assert_eq!(physical.before.first_mass_moment_kg_m, [15.0, 22.5, 30.0]);
    assert_eq!(physical.before.momentum_kg_m_s, [2.0, -3.0, 4.0]);
    assert_eq!(
        physical.before.angular_momentum_kg_m2_s,
        [-12.0, 15.0, 18.0]
    );
    assert_eq!(physical.before.kinetic_energy_j, 0.07 * 40.0);
    assert_eq!(physical.kinetic_energy_change_j, -6e-15 * 40.0);
    assert_eq!(
        physical.absolute_angular_momentum_defect_kg_m2_s,
        5e-15 * 30.0
    );
    assert_eq!(physical.maximum_viscous_stress_defect_pa, 7e-15 * 16.0);
    assert_eq!(physical.before, physical.after);
    assert_eq!(physical.source_geometry_epoch, 9);
    assert_eq!(physical.target_geometry_epoch, 10);
}

#[test]
fn rounded_si_endpoint_inventory_cannot_hide_behind_a_small_native_receipt() {
    let mut native = report();
    native.absolute_mass_defect = 0.0;
    native.before.mass = 1e12;
    native.after.mass = native.before.mass.next_up();
    let mut budgets = [1.0; 6];
    budgets[0] = 0.0;
    assert!(receipt(native, [1.0; 6], budgets).is_err());
    native.after = native.before;
    native.maximum_viscous_stress_defect = 2.0;
    assert!(receipt(native, [1.0; 6], [1.0; 6]).is_err());
}

#[test]
fn si_underflow_overflow_and_unrepresentable_angular_inventory_refuse() {
    assert!(times(f64::from_bits(1), 0.5).is_err());
    assert!(times(f64::MAX, 2.0).is_err());
    assert_eq!(times(0.0, 2.0).unwrap(), 0.0);
    let mut native = report();
    native.before.angular_momentum[0] = f64::MAX;
    assert!(receipt(native, [2.0; 6], [f64::MAX; 6]).is_err());
    native = report();
    native.kinetic_energy_change = f64::from_bits(1);
    assert!(receipt(native, [0.5; 6], [1.0; 6]).is_err());
}
