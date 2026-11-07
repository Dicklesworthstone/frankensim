//! Physical SDIRK2 stages on the native full-tensor CHT face operator.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance, SolidMaterial,
    TENSOR_SDIRK2_GAMMA, ThermalFace, ThermalSetup, Voxel, VoxelDomain, advance_tensor_conduction,
    advance_tensor_conduction_sdirk2,
};

fn diagonal_tensor() -> ConductivityTensor {
    ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]).unwrap()
}

fn rotated_tensor() -> ConductivityTensor {
    ConductivityTensor::new([[4.0, 1.0, 0.4], [1.0, 2.0, 0.2], [0.4, 0.2, 3.0]]).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-8, "{actual} != {expected}");
}

fn pair_setup(resistance_m2_k_w: f64) -> ThermalSetup {
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w,
    });
    setup
}

#[test]
fn insulated_finite_contact_pair_retains_both_stage_heats_and_original_storage() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(); 2];
    let capacity = [2.0, 3.0];
    let old = [400.0, 300.0];
    let setup = pair_setup(0.5);
    let dt = 2.0;
    let step = advance_tensor_conduction_sdirk2(
        &domain,
        &solids,
        [&setup, &setup],
        &tensors,
        &capacity,
        &old,
        5.0,
        dt,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    let gamma = TENSOR_SDIRK2_GAMMA;
    let a = gamma * dt;
    // Unit area, two half-cell bulk resistances, and one finite contact.
    let g = 1.0 / (0.5 + 0.5 + 0.5);
    let lambda = g * (1.0 / capacity[0] + 1.0 / capacity[1]);
    let differences = [
        100.0 / (1.0 + a * lambda),
        100.0 * (1.0 + (2.0 * gamma - 1.0) * lambda * dt) / (1.0 + a * lambda).powi(2),
    ];
    let original_energy = capacity[0] * old[0] + capacity[1] * old[1];
    for (stage, difference) in step.stages.iter().zip(differences) {
        close(
            stage.temperature[0],
            (original_energy + capacity[1] * difference) / 5.0,
        );
        close(
            stage.temperature[1],
            (original_energy - capacity[0] * difference) / 5.0,
        );
        let heat = g * difference;
        close(stage.outward_heat_w[0][1], heat);
        close(stage.outward_heat_w[1][0], -heat);
        close(
            stage.face_temperature_k[0][1] - stage.face_temperature_k[1][0],
            0.5 * heat,
        );
        assert_eq!(stage.report.balance.source_w, 0.0);
    }
    close(step.record.stored_energy_change_j, 0.0);
    close(step.record.source_j, 0.0);
    close(step.record.boundary_outflow_j, 0.0);
    close(step.record.closure_j, 0.0);
    assert_eq!(step.stage_time_s, [5.0 + a, 7.0]);
    assert_eq!(step.record.time_s, 7.0);
    assert_eq!(
        step.record.iterations,
        step.stages[0].report.iterations + step.stages[1].report.iterations
    );
    assert_eq!(old, [400.0, 300.0]);
}

#[test]
fn native_finite_contact_decay_converges_at_second_order() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(); 2];
    let capacity = [2.0, 3.0];
    let setup = pair_setup(0.5);
    let config = EnergyConfig::default();
    let gate = CancelGate::new();
    let march = |steps: usize| {
        let dt = 1.0 / steps as f64;
        let mut previous = vec![400.0, 300.0];
        for index in 0..steps {
            let accepted = advance_tensor_conduction_sdirk2(
                &domain,
                &solids,
                [&setup, &setup],
                &tensors,
                &capacity,
                &previous,
                index as f64 * dt,
                dt,
                &config,
                &gate,
            )
            .unwrap();
            close(accepted.record.closure_j, 0.0);
            let [_, endpoint] = accepted.stages;
            previous = endpoint.temperature;
        }
        previous[0] - previous[1]
    };
    let lambda: f64 = (2.0 / 3.0) * (1.0 / 2.0 + 1.0 / 3.0);
    let exact = 100.0 * (-lambda).exp();
    let coarse = (march(10) - exact).abs();
    let fine = (march(20) - exact).abs();
    assert!((3.9..4.1).contains(&(coarse / fine)), "{coarse} {fine}");
}

#[test]
fn rotated_tensor_linear_forcing_uses_each_tableau_sample_and_capacity_once() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 0.25, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("rotated", 1.0)];
    let tensors = [rotated_tensor(); 3];
    let capacity = [1.0, 2.0, 4.0]; // Complete J/K despite non-unit cell volume.
    let start = 3.0;
    let dt = 2.0;
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let rate = |t: f64| 5.0 + 7.0 * t;
    let mut first = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let mut second = first.clone();
    first.power_w = capacity.map(|c| c * rate(start + a)).to_vec();
    second.power_w = capacity.map(|c| c * rate(start + dt)).to_vec();
    let step = advance_tensor_conduction_sdirk2(
        &domain,
        &solids,
        [&first, &second],
        &tensors,
        &capacity,
        &[300.0; 3],
        start,
        dt,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    let change = 5.0 * dt + 3.5 * ((start + dt).powi(2) - start.powi(2));
    for cell in 0..3 {
        close(
            step.stages[0].temperature[cell],
            300.0 + a * rate(start + a),
        );
        close(step.stages[1].temperature[cell], 300.0 + change);
        for stage in &step.stages {
            for &heat in &stage.outward_heat_w[cell] {
                close(heat, 0.0);
            }
        }
    }
    // The internal carry must not become generated power in the native ledger.
    close(
        step.stages[0].report.balance.source_w,
        7.0 * rate(start + a),
    );
    close(
        step.stages[1].report.balance.source_w,
        7.0 * rate(start + dt),
    );
    close(step.record.source_j, 7.0 * change);
    close(step.record.stored_energy_change_j, 7.0 * change);
    close(step.record.boundary_outflow_j, 0.0);
    close(step.record.closure_j, 0.0);
}

#[test]
fn changed_fixed_masks_and_contacts_integrate_each_physical_reservoir_once() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(); 2];
    let capacity = [2.0, 3.0];
    let old = [300.0, 310.0];
    let mut first = pair_setup(0.5);
    first.fixed_temperature.push((0, 350.0));
    first.power_w = vec![1000.0, 2.0];
    let mut second = pair_setup(0.25);
    second.fixed_temperature.push((1, 330.0));
    second.power_w = vec![3.0, 2000.0];
    let dt = 2.0;
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
    let r = b / a;
    let g1 = 1.0 / 1.5;
    let g2 = 1.0 / 1.25;
    let t1 = [
        350.0,
        (capacity[1] * old[1] + a * (2.0 + g1 * 350.0)) / (capacity[1] + a * g1),
    ];
    let q1 = g1 * (t1[0] - t1[1]);
    let rate1 = [capacity[0] * (t1[0] - old[0]) / a, 2.0 + q1];
    let t2 = [
        (capacity[0] * old[0] + a * (3.0 + r * rate1[0] + g2 * 330.0)) / (capacity[0] + a * g2),
        330.0,
    ];
    let q2 = g2 * (t2[0] - t2[1]);
    let reservoir1 = q1 + rate1[0];
    let reservoir2 = -q2 + (capacity[1] * (t2[1] - old[1]) - b * rate1[1]) / a;
    let step = advance_tensor_conduction_sdirk2(
        &domain,
        &solids,
        [&first, &second],
        &tensors,
        &capacity,
        &old,
        0.0,
        dt,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    for cell in 0..2 {
        close(step.stages[0].temperature[cell], t1[cell]);
        close(step.stages[1].temperature[cell], t2[cell]);
    }
    close(step.stages[0].outward_heat_w[0][1], q1);
    close(step.stages[1].outward_heat_w[0][1], q2);
    close(
        step.stages[0].face_temperature_k[0][1] - step.stages[0].face_temperature_k[1][0],
        0.5 * q1,
    );
    close(
        step.stages[1].face_temperature_k[0][1] - step.stages[1].face_temperature_k[1][0],
        0.25 * q2,
    );
    close(step.record.source_j, b * 2.0 + a * 3.0);
    close(
        step.record.boundary_outflow_j,
        -b * reservoir1 - a * reservoir2,
    );
    close(
        step.record.stored_energy_change_j,
        capacity[0] * (t2[0] - old[0]) + capacity[1] * (t2[1] - old[1]),
    );
    close(step.record.closure_j, 0.0);
    assert_eq!(old, [300.0, 310.0]);
    assert_eq!(first.fixed_temperature, vec![(0, 350.0)]);
    assert_eq!(second.fixed_temperature, vec![(1, 330.0)]);
}

#[test]
fn prescribed_body_endpoint_storage_does_not_double_count_stage_one() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let mut first = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let mut second = first.clone();
    first.fixed_temperature.push((0, 340.0));
    second.fixed_temperature.push((0, 360.0));
    let step = advance_tensor_conduction_sdirk2(
        &domain,
        &solids,
        [&first, &second],
        &[rotated_tensor()],
        &[7.0],
        &[300.0],
        1.0,
        3.0,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    close(step.stages[0].temperature[0], 340.0);
    close(step.stages[1].temperature[0], 360.0);
    close(step.record.stored_energy_change_j, 420.0);
    close(step.record.source_j, 0.0);
    close(step.record.boundary_outflow_j, -420.0);
    close(step.record.closure_j, 0.0);
}

#[test]
fn krylov_work_is_shared_and_failures_publish_no_partial_interval() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("body", 1.0)];
    let tensors = [diagonal_tensor(); 2];
    let capacity = [2.0, 3.0];
    let old = [400.0, 300.0];
    let setup = pair_setup(0.5);
    let config = EnergyConfig::default();
    let first = advance_tensor_conduction(
        &domain,
        &solids,
        &setup,
        &tensors,
        &capacity,
        &old,
        0.0,
        TENSOR_SDIRK2_GAMMA,
        &config,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(first.solution.report.iterations > 0);
    let limited = EnergyConfig {
        max_iterations: first.solution.report.iterations,
        ..config
    };
    let solve = |second: &ThermalSetup, config: &EnergyConfig, gate: &CancelGate| {
        advance_tensor_conduction_sdirk2(
            &domain,
            &solids,
            [&setup, second],
            &tensors,
            &capacity,
            &old,
            0.0,
            1.0,
            config,
            gate,
        )
    };
    assert!(matches!(
        solve(&setup, &limited, &CancelGate::new()),
        Err(ChtError::SolverNotConverged { .. })
    ));
    let mut invalid_endpoint = setup.clone();
    invalid_endpoint
        .fixed_temperature
        .push((domain.cell_count(), 350.0));
    assert!(matches!(
        solve(&invalid_endpoint, &config, &CancelGate::new()),
        Err(ChtError::InvalidInput { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        solve(&setup, &config, &gate).unwrap_err(),
        ChtError::Cancelled
    );
    assert_eq!(old, [400.0, 300.0]);
    assert!(setup.fixed_temperature.is_empty());
    assert_eq!(setup.face_contacts.len(), 1);
}
