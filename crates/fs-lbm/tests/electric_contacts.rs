//! Synthetic circuits through the production voxel current and heat owners.
//! These are equation regressions, not identified or measured contact data.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, Conductor, ElectricalContact, ElectricSetup, Electrode, ElectrodeDrive,
    Voxel, VoxelDomain, solve_electric,
};

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1.0), "{a} != {b}");
}
fn bar(drive: ElectrodeDrive) -> (VoxelDomain, ElectricSetup) {
    let domain = VoxelDomain::from_fn(5, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let setup = ElectricSetup {
        conductors: vec![Some(Conductor::new(2.0))],
        electrodes: vec![
            Electrode { lo: [0, 0, 0], hi: [1, 1, 1], drive: ElectrodeDrive::Voltage(0.0) },
            Electrode { lo: [4, 0, 0], hi: [5, 1, 1], drive },
        ],
    };
    (domain, setup)
}
fn contacts() -> [ElectricalContact; 2] {
    [ElectricalContact { cells: [1, 2], resistance_ohm_m2: 4.0, heat_fraction_to_first: 0.25 },
     ElectricalContact { cells: [3, 4], resistance_ohm_m2: 2.0, heat_fraction_to_first: 0.25 }]
}

#[test]
fn bulk_and_contact_resistance_set_current_and_each_cells_heat() {
    for drive in [ElectrodeDrive::Voltage(9.0), ElectrodeDrive::Current(0.75),
        ElectrodeDrive::Current(-0.75)]
    {
        let (domain, setup) = bar(drive);
        let field = setup.solve_with_contacts(
            &domain, &[300.0; 5], &contacts(), |_, _, _| Ok(2.0), &CancelGate::new(),
        ).unwrap();
        let sign = if matches!(drive, ElectrodeDrive::Current(i) if i < 0.0) { -1.0 } else { 1.0 };
        // Bulk resistance 3*2 plus contacts 4+2 = 12 ohm.
        near(field.electrodes[1].voltage_v, sign * 9.0);
        near(field.electrodes[1].current_a, sign * 0.75);
        near(field.electrodes[0].current_a, -sign * 0.75);
        for (got, expected) in field.joule_w.iter().zip([0.0, 1.6875, 2.8125, 1.40625, 0.84375]) {
            near(*got, expected);
        }
        near(field.total_joule_w, 6.75);
        near(field.delivered_w, field.total_joule_w);
    }
}

#[test]
fn heat_partition_changes_only_interface_deposition_not_electrical_response() {
    let (domain, setup) = bar(ElectrodeDrive::Current(1.0));
    let solve = |fraction| setup.solve_with_contacts(
        &domain, &[300.0; 5], &[ElectricalContact {
            cells: [1, 2], resistance_ohm_m2: 4.0, heat_fraction_to_first: fraction,
        }], |_, _, _| Ok(2.0), &CancelGate::new(),
    ).unwrap();
    let a = solve(0.0);
    let b = solve(1.0);
    assert_eq!(a.potential_v, b.potential_v);
    assert_eq!(a.electrodes, b.electrodes);
    near(a.total_joule_w, 10.0);
    near(b.total_joule_w, a.total_joule_w);
    near(a.joule_w[1], 2.0);
    near(a.joule_w[2], 6.0);
    near(b.joule_w[1], 6.0);
    near(b.joule_w[2], 2.0);
}

#[test]
fn empty_and_zero_contacts_reproduce_the_original_solver_bitwise() {
    let (domain, setup) = bar(ElectrodeDrive::Voltage(9.0));
    let gate = CancelGate::new();
    let original = solve_electric(&domain, &setup, Some(&[300.0; 5]), &gate).unwrap();
    for input in [vec![], contacts().map(|mut c| { c.resistance_ohm_m2 = 0.0; c }).to_vec()] {
        let result = setup.solve_with_contacts(
            &domain, &[300.0; 5], &input, |_, _, _| Ok(2.0), &gate,
        ).unwrap();
        assert_eq!(result, original);
    }
}

#[test]
fn reversing_surface_order_and_input_order_preserves_the_physical_result() {
    let (domain, setup) = bar(ElectrodeDrive::Voltage(9.0));
    let original = contacts();
    let mut reversed = original.map(|face| ElectricalContact {
        cells: [face.cells[1], face.cells[0]],
        heat_fraction_to_first: 1.0 - face.heat_fraction_to_first, ..face
    });
    reversed.reverse();
    let run = |faces: &[ElectricalContact]| setup.solve_with_contacts(
        &domain, &[300.0; 5], faces, |_, _, _| Ok(2.0), &CancelGate::new(),
    ).unwrap();
    assert_eq!(run(&original), run(&reversed));
}

#[test]
fn contact_specific_resistance_uses_actual_face_area_on_a_graded_grid() {
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.4, 1.0, 2.0, 0.5], vec![2.0], vec![3.0]], |_| Voxel::Solid(0),
    ).unwrap();
    let setup = ElectricSetup {
        conductors: vec![Some(Conductor::new(1.0))],
        electrodes: vec![
            Electrode { lo: [0, 0, 0], hi: [1, 1, 1], drive: ElectrodeDrive::Voltage(0.0) },
            Electrode { lo: [3, 0, 0], hi: [4, 1, 1], drive: ElectrodeDrive::Voltage(10.0) },
        ],
    };
    let field = setup.solve_with_contacts(
        &domain, &[300.0; 4], &[ElectricalContact {
            cells: [1, 2], resistance_ohm_m2: 12.0, heat_fraction_to_first: 0.25,
        }], |cell, _, _| Ok(if cell == 1 { 2.0 } else { 4.0 }), &CancelGate::new(),
    ).unwrap();
    let current = 10.0 / ((2.0 * 1.0 + 4.0 * 2.0 + 12.0) / 6.0);
    near(field.electrodes[1].current_a, current);
    near(field.joule_w[1], current * current * (2.0 + 0.25 * 12.0) / 6.0);
    near(field.joule_w[2], current * current * (8.0 + 0.75 * 12.0) / 6.0);
    near(field.total_joule_w, 10.0 * current);
}

#[test]
fn one_contact_redistributes_current_across_parallel_paths() {
    let domain = VoxelDomain::from_fn(3, 2, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let setup = ElectricSetup {
        conductors: vec![Some(Conductor::new(1.0))],
        electrodes: vec![
            Electrode { lo: [0, 0, 0], hi: [1, 2, 1], drive: ElectrodeDrive::Voltage(0.0) },
            Electrode { lo: [2, 0, 0], hi: [3, 2, 1], drive: ElectrodeDrive::Voltage(1.0) },
        ],
    };
    let field = setup.solve_with_contacts(
        &domain, &[300.0; 6], &[ElectricalContact {
            cells: [1, 2], resistance_ohm_m2: 1.5, heat_fraction_to_first: 0.0,
        }], |_, _, _| Ok(1.0), &CancelGate::new(),
    ).unwrap();
    // The two free rows solve [[3.5,-1],[-1,5]] v = [0.5,2].
    near(field.potential_v[1], 3.0 / 11.0);
    near(field.potential_v[4], 5.0 / 11.0);
    near(field.electrodes[1].current_a, 16.0 / 11.0);
    near(field.joule_w[2], 24.0 / 121.0); // contact heat in the ideal electrode
    near(field.total_joule_w, field.delivered_w);
}

#[test]
fn invalid_contacts_and_late_cancellation_publish_nothing_and_preserve_inputs() {
    let (mut domain, setup) = bar(ElectrodeDrive::Voltage(9.0));
    let before = setup.clone();
    let first = contacts()[0];
    let mut invalid = vec![vec![first, first], vec![first, ElectricalContact {
        cells: [2, 1], ..first
    }]];
    for cells in [[1, 1], [0, 4], [1, usize::MAX]] {
        invalid.push(vec![ElectricalContact { cells, ..first }]);
    }
    for resistance in [-1.0, f64::NAN, f64::INFINITY] {
        invalid.push(vec![ElectricalContact { resistance_ohm_m2: resistance, ..first }]);
    }
    for share in [-0.1, 1.1, f64::NAN] {
        invalid.push(vec![ElectricalContact { heat_fraction_to_first: share, ..first }]);
    }
    for faces in invalid {
        let mut calls = 0;
        assert!(matches!(setup.solve_with_contacts(
            &domain, &[300.0; 5], &faces, |_, _, _| { calls += 1; Ok(2.0) }, &CancelGate::new(),
        ), Err(ChtError::InvalidInput { field: "electric.contacts", .. })));
        assert!(calls <= 3);
    }
    let gate = CancelGate::new();
    assert_eq!(setup.solve_with_contacts(
        &domain, &[300.0; 5], &contacts(), |cell, _, _| {
            if cell == 3 { gate.request(); }
            Ok(2.0)
        }, &gate,
    ), Err(ChtError::Cancelled));
    domain.set(2, 0, 0, Voxel::Fluid);
    assert!(matches!(setup.solve_with_contacts(
        &domain, &[300.0; 5], &[first], |_, _, _| Ok(2.0), &CancelGate::new(),
    ), Err(ChtError::InvalidInput { field: "electric.contacts", .. })));
    assert_eq!(setup, before);
}
