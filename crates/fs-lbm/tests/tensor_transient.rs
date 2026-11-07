use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ConductivityTensor, EnergyConfig, FaceContactResistance, SolidMaterial, ThermalFace,
    ThermalSetup, Voxel, VoxelDomain, advance_tensor_conduction,
};

fn diagonal_tensor() -> ConductivityTensor {
    ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]).unwrap()
}

fn rotated_tensor() -> ConductivityTensor {
    ConductivityTensor::new([[4.0, 1.0, 0.4], [1.0, 2.0, 0.2], [0.4, 0.2, 3.0]]).unwrap()
}

#[test]
fn insulated_tensor_pair_uses_cell_capacities_and_conserves_contact_heat() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(), diagonal_tensor()];
    let capacity = [2.0, 3.0];
    let previous = [400.0, 300.0];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.5,
    });
    let dt = 2.0;
    let next = advance_tensor_conduction(
        &domain,
        &solids,
        &setup,
        &tensors,
        &capacity,
        &previous,
        5.0,
        dt,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    // Two half-cell bulk resistances plus the original finite joint.
    let g = 1.0 / (0.5 + 0.5 + 0.5);
    let difference = 100.0 / (1.0 + dt * g * (1.0 / capacity[0] + 1.0 / capacity[1]));
    let energy = capacity[0] * previous[0] + capacity[1] * previous[1];
    let expected = [
        (energy + capacity[1] * difference) / (capacity[0] + capacity[1]),
        (energy - capacity[0] * difference) / (capacity[0] + capacity[1]),
    ];
    for (&actual, &reference) in next.solution.temperature.iter().zip(&expected) {
        assert!((actual - reference).abs() < 2e-8, "{actual} != {reference}");
    }
    let heat = g * difference;
    assert!((next.solution.outward_heat_w[0][1] - heat).abs() < 2e-8);
    assert!((next.solution.outward_heat_w[1][0] + heat).abs() < 2e-8);
    let jump = next.solution.face_temperature_k[0][1] - next.solution.face_temperature_k[1][0];
    assert!((jump - 0.5 * heat).abs() < 2e-8);
    assert!(next.record.stored_energy_change_j.abs() < 2e-8);
    assert!(next.record.closure_j.abs() < 2e-8);
    assert_eq!(next.record.time_s, 7.0);
    assert_eq!(previous, [400.0, 300.0]);
}

#[test]
fn rotated_tensor_uniform_heating_stores_each_capacity_once() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 0.25, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("anisotropic body", 1.0)];
    let tensors = [rotated_tensor(), rotated_tensor(), rotated_tensor()];
    // Explicit J/K, despite the non-unit voxel volume.
    let capacity = [1.0, 2.0, 4.0];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![5.0, 10.0, 20.0];
    let next = advance_tensor_conduction(
        &domain,
        &solids,
        &setup,
        &tensors,
        &capacity,
        &[300.0; 3],
        0.0,
        2.0,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    for &temperature in &next.solution.temperature {
        assert!((temperature - 310.0).abs() < 2e-8);
    }
    for faces in &next.solution.outward_heat_w {
        assert!(faces.iter().all(|q| q.abs() < 2e-8));
    }
    assert!((next.record.stored_energy_change_j - 70.0).abs() < 2e-8);
    assert!((next.record.source_j - 70.0).abs() < 2e-8);
    assert!(next.record.boundary_outflow_j.abs() < 2e-8);
    assert!(next.record.closure_j.abs() < 2e-8);
}

#[test]
fn fixed_tensor_cell_reservoir_supplies_the_storage_jump_once() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature.push((0, 350.0));
    let next = advance_tensor_conduction(
        &domain,
        &solids,
        &setup,
        &[rotated_tensor()],
        &[7.0],
        &[300.0],
        8.0,
        2.0,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!((next.solution.temperature[0] - 350.0).abs() < 2e-8);
    assert!((next.record.stored_energy_change_j - 350.0).abs() < 2e-8);
    assert!((next.record.boundary_outflow_j + 350.0).abs() < 2e-8);
    assert_eq!(next.record.source_j, 0.0);
    assert!(next.record.closure_j.abs() < 2e-8);
    assert_eq!(next.record.time_s, 10.0);
}

#[test]
fn tensor_pair_backward_euler_converges_to_the_physical_decay_mode() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(), diagonal_tensor()];
    let capacity = [2.0, 3.0];
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let config = EnergyConfig::default();
    let gate = CancelGate::new();
    let march = |steps: usize| {
        let dt = 1.0 / steps as f64;
        let mut temperature = vec![400.0, 300.0];
        for index in 0..steps {
            let next = advance_tensor_conduction(
                &domain,
                &solids,
                &setup,
                &tensors,
                &capacity,
                &temperature,
                index as f64 * dt,
                dt,
                &config,
                &gate,
            )
            .unwrap();
            assert!(next.record.closure_j.abs() < 2e-8);
            temperature = next.solution.temperature;
        }
        temperature[0] - temperature[1]
    };
    // Unit area, unit total bulk resistance, and no extra joint.
    let exact = 100.0 * (-(1.0_f64 / 2.0 + 1.0 / 3.0)).exp();
    let coarse = (march(10) - exact).abs();
    let fine = (march(20) - exact).abs();
    assert!(
        coarse / fine > 1.9 && coarse / fine < 2.1,
        "{coarse} {fine}"
    );
}
