//! Original conductivity integrals and physical SDIRK2 heat/storage accounting.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CellSink, ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance,
    RotatedSolidConductivity, SolidMaterial, TENSOR_SDIRK2_GAMMA, TensorConductivityLaw,
    TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_nonlinear_tensor_conduction_sdirk2, advance_tensor_conduction_sdirk2,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        ..TensorNonlinearConfig::default()
    }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-7, "{actual} != {expected}");
}

fn linear_material() -> SolidMaterial {
    SolidMaterial::new("original k(T)=T/300", 1.0)
        .with_conductivity_table(&[(200.0, 2.0 / 3.0), (600.0, 2.0)])
}

fn cooling_setup() -> ThermalSetup {
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    ThermalSetup::new(faces)
}

#[test]
fn nonlinear_stages_resolve_exact_integrals_and_original_state_storage() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_material()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let setup = cooling_setup();
    let previous = [400.0];
    let capacity = [2.0];
    let dt = 0.8;
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
    let step = advance_nonlinear_tensor_conduction_sdirk2(
        &domain,
        &material,
        [&setup, &setup],
        &[&law],
        &capacity,
        &previous,
        4.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // The physical half-cell law is Q=(T²-300²)/300, not k(T)*(T-300)/.5.
    // Each stage solves C*(T-T0)/a + Q = carry. Solve its quadratic directly.
    let temperature = |carry: f64| {
        let linear = 300.0 * capacity[0] / a;
        let constant = linear * previous[0] + 300.0_f64.powi(2) + 300.0 * carry;
        2.0 * constant / (linear + (linear * linear + 4.0 * constant).sqrt())
    };
    let heat = |t: f64| (t * t - 300.0_f64.powi(2)) / 300.0;
    let first = temperature(0.0);
    let second = temperature(-b / a * heat(first));
    for (stage, expected) in step.stages.iter().zip([first, second]) {
        close(stage.solution.temperature[0], expected);
        close(stage.solution.face_temperature_k[0][0], 300.0);
        close(stage.solution.outward_heat_w[0][0], heat(expected));
        for face in 1..6 {
            close(stage.solution.outward_heat_w[0][face], 0.0);
        }
        close(stage.solution.report.balance.source_w, 0.0);
        assert!(stage.newton_iterations > 0);
    }
    close(
        step.record.stored_energy_change_j,
        capacity[0] * (second - previous[0]),
    );
    close(
        step.record.boundary_outflow_j,
        b * heat(first) + a * heat(second),
    );
    close(step.record.source_j, 0.0);
    close(step.record.closure_j, 0.0);
    assert_eq!(step.stage_time_s, [4.0 + a, 4.0 + dt]);
    assert_eq!(step.record.time_s, 4.0 + dt);
    assert_eq!(
        step.record.iterations,
        step.stages
            .iter()
            .map(|s| s.solution.report.iterations)
            .sum::<usize>()
    );
    assert_eq!(previous, [400.0]);
}

#[test]
fn original_nonlinear_cooling_converges_at_second_order_in_time() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_material()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let setup = cooling_setup();
    let config = controls();
    let gate = CancelGate::new();
    let march = |count: usize| {
        let dt = 1.0 / count as f64;
        let mut previous = vec![400.0];
        for index in 0..count {
            let step = advance_nonlinear_tensor_conduction_sdirk2(
                &domain,
                &material,
                [&setup, &setup],
                &[&law],
                &[2.0],
                &previous,
                index as f64 * dt,
                dt,
                &config,
                &gate,
            )
            .unwrap();
            close(step.record.closure_j, 0.0);
            let [_, endpoint] = step.stages;
            previous = endpoint.solution.temperature;
        }
        previous[0]
    };
    // C*dT/dt=-(T²-Tb²)/300, so (T-Tb)/(T+Tb) decays as
    // exp(-2*Tb*t/(300*C)). This is the continuous nonlinear solution.
    let ratio = (100.0 / 700.0) * (-1.0_f64).exp();
    let exact = 300.0 * (1.0 + ratio) / (1.0 - ratio);
    let coarse = (march(20) - exact).abs();
    let fine = (march(40) - exact).abs();
    assert!((3.9..4.15).contains(&(coarse / fine)), "{coarse} {fine}");
    assert!(fine < 0.01, "fine-step temperature error: {fine}");
}

#[test]
fn nonlinear_finite_contact_preserves_each_stage_law_and_total_energy() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_material()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.5,
    });
    let previous = [400.0, 300.0];
    let capacity = [2.0, 3.0];
    let dt = 2.0;
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let r = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
    let step = advance_nonlinear_tensor_conduction_sdirk2(
        &domain,
        &material,
        [&setup, &setup],
        &[&law, &law],
        &capacity,
        &previous,
        0.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut carried = [0.0; 2];
    for stage in &step.stages {
        let solved = &stage.solution;
        let t = &solved.temperature;
        let surface = [
            solved.face_temperature_k[0][1],
            solved.face_temperature_k[1][0],
        ];
        let q = solved.outward_heat_w[0][1];
        assert!(q > 0.0);
        close(solved.outward_heat_w[1][0], -q);
        close((t[0] * t[0] - surface[0] * surface[0]) / 300.0, q);
        close((surface[1] * surface[1] - t[1] * t[1]) / 300.0, q);
        close(surface[0] - surface[1], 0.5 * q);
        close(capacity[0] * (t[0] - previous[0]) / a, carried[0] - q);
        close(capacity[1] * (t[1] - previous[1]) / a, carried[1] + q);
        carried = [-r * q, r * q];
    }
    close(step.record.stored_energy_change_j, 0.0);
    close(step.record.source_j, 0.0);
    close(step.record.boundary_outflow_j, 0.0);
    close(step.record.closure_j, 0.0);
    assert_eq!(previous, [400.0, 300.0]);
}

#[test]
fn constant_rotated_laws_match_frozen_stages_with_graded_geometry_and_films() {
    let domain = VoxelDomain::graded_from_fn([vec![0.1, 0.3, 0.6], vec![0.2], vec![0.3]], |_| {
        Voxel::Solid(0)
    })
    .unwrap();
    let material = [SolidMaterial::new("constant owner", 1.0)];
    let tensor =
        ConductivityTensor::new([[5.0, 1.0, 0.3], [1.0, 3.0, 0.2], [0.3, 0.2, 2.0]]).unwrap();
    let mut first = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    first.faces[0] = ThermalFace::Temperature(400.0);
    first.faces[1] = ThermalFace::Convective {
        h: 8.0,
        ambient: 330.0,
    };
    first.power_w = vec![1.0, 0.5, 2.0];
    first.cell_sinks.push(CellSink {
        cell: 1,
        conductance_w_k: 0.4,
        temperature_k: 310.0,
    });
    first.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.05,
    });
    let mut second = first.clone();
    second.faces[0] = ThermalFace::Temperature(390.0);
    second.faces[1] = ThermalFace::Convective {
        h: 8.0,
        ambient: 325.0,
    };
    second.power_w = vec![2.0, 0.75, 1.5];
    second.cell_sinks[0].temperature_k = 315.0;
    second.face_contacts[0].resistance_m2_k_w = 0.08;
    let capacity = [0.5, 2.0, 3.0];
    let previous = [350.0, 330.0, 320.0];
    let config = controls();
    let gate = CancelGate::new();
    let frozen = advance_tensor_conduction_sdirk2(
        &domain,
        &material,
        [&first, &second],
        &[tensor; 3],
        &capacity,
        &previous,
        1.0,
        0.4,
        &config.linear,
        &gate,
    )
    .unwrap();
    let nonlinear = advance_nonlinear_tensor_conduction_sdirk2(
        &domain,
        &material,
        [&first, &second],
        &[&tensor, &tensor, &tensor],
        &capacity,
        &previous,
        1.0,
        0.4,
        &config,
        &gate,
    )
    .unwrap();
    for (actual, expected) in nonlinear.stages.iter().zip(&frozen.stages) {
        for cell in 0..3 {
            close(
                actual.solution.temperature[cell],
                expected.temperature[cell],
            );
            for face in 0..6 {
                close(
                    actual.solution.face_temperature_k[cell][face],
                    expected.face_temperature_k[cell][face],
                );
                close(
                    actual.solution.outward_heat_w[cell][face],
                    expected.outward_heat_w[cell][face],
                );
            }
        }
    }
    close(
        nonlinear.record.stored_energy_change_j,
        frozen.record.stored_energy_change_j,
    );
    close(nonlinear.record.source_j, frozen.record.source_j);
    close(
        nonlinear.record.boundary_outflow_j,
        frozen.record.boundary_outflow_j,
    );
    close(nonlinear.record.closure_j, 0.0);
}

#[test]
fn changed_prescribed_masks_count_original_nonlinear_reservoir_heat_once() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_material()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let mut first = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    first.fixed_temperature.push((0, 350.0));
    first.power_w = vec![1000.0, 2.0];
    first.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.5,
    });
    let mut second = first.clone();
    second.fixed_temperature = vec![(1, 330.0)];
    second.power_w = vec![3.0, 2000.0];
    second.face_contacts[0].resistance_m2_k_w = 0.25;
    let capacity = [2.0, 3.0];
    let previous = [300.0, 310.0];
    let dt = 2.0;
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
    let step = advance_nonlinear_tensor_conduction_sdirk2(
        &domain,
        &material,
        [&first, &second],
        &[&law, &law],
        &capacity,
        &previous,
        0.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let t1 = &step.stages[0].solution.temperature;
    let t2 = &step.stages[1].solution.temperature;
    let q1 = step.stages[0].solution.outward_heat_w[0][1];
    let q2 = step.stages[1].solution.outward_heat_w[0][1];
    close(t1[0], 350.0);
    close(t2[1], 330.0);
    close(capacity[1] * (t1[1] - previous[1]) / a, 2.0 + q1);
    let carry0 = b / a * capacity[0] * (t1[0] - previous[0]) / a;
    let carry1 = b / a * (2.0 + q1);
    close(capacity[0] * (t2[0] - previous[0]) / a, 3.0 - q2 + carry0);
    let reservoir1 = -q1 - capacity[0] * (t1[0] - previous[0]) / a;
    let reservoir2 = q2 - capacity[1] * (t2[1] - previous[1]) / a + carry1;
    close(step.record.source_j, b * 2.0 + a * 3.0);
    close(
        step.record.boundary_outflow_j,
        b * reservoir1 + a * reservoir2,
    );
    close(
        step.record.stored_energy_change_j,
        capacity[0] * (t2[0] - previous[0]) + capacity[1] * (t2[1] - previous[1]),
    );
    close(step.record.closure_j, 0.0);
    assert_eq!(previous, [300.0, 310.0]);
    assert_eq!(first.fixed_temperature, vec![(0, 350.0)]);
    assert_eq!(second.fixed_temperature, vec![(1, 330.0)]);
}

#[test]
fn shared_nonlinear_and_linear_budgets_and_failed_second_stage_publish_nothing() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_material()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let setup = cooling_setup();
    let previous = [400.0];
    let config = controls();
    let solve = |second: &ThermalSetup, config: &TensorNonlinearConfig, gate: &CancelGate| {
        advance_nonlinear_tensor_conduction_sdirk2(
            &domain,
            &material,
            [&setup, second],
            &[&law],
            &[2.0],
            &previous,
            0.0,
            1.0,
            config,
            gate,
        )
    };
    let accepted = solve(&setup, &config, &CancelGate::new()).unwrap();
    let newton: usize = accepted.stages.iter().map(|s| s.newton_iterations).sum();
    assert!(accepted.stages.iter().all(|s| s.newton_iterations > 0));
    assert!(
        accepted
            .stages
            .iter()
            .all(|s| s.solution.report.iterations > 0)
    );
    for limited in [
        TensorNonlinearConfig {
            max_newton_iterations: newton - 1,
            ..config
        },
        TensorNonlinearConfig {
            linear: EnergyConfig {
                max_iterations: accepted.record.iterations - 1,
                ..config.linear
            },
            ..config
        },
    ] {
        assert!(matches!(
            solve(&setup, &limited, &CancelGate::new()),
            Err(ChtError::SolverNotConverged { .. })
        ));
    }
    let mut invalid_endpoint = setup.clone();
    invalid_endpoint.fixed_temperature.push((1, 350.0));
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
    assert_eq!(previous, [400.0]);
    assert_eq!(setup, cooling_setup());
}

struct BoundedOriginalLaw(ConductivityTensor);

impl TensorConductivityLaw for BoundedOriginalLaw {
    fn conductivity_at_temperature(&self, t: f64) -> Result<ConductivityTensor, ChtError> {
        if t.is_finite() && (250.0..=450.0).contains(&t) {
            Ok(self.0)
        } else {
            Err(ChtError::InvalidInput {
                field: "test.bounded_original_law",
                reason: "original conductivity support is 250..=450 K".into(),
            })
        }
    }

    fn mean_conductivity(&self, from: f64, to: f64) -> Result<ConductivityTensor, ChtError> {
        self.conductivity_at_temperature(from)?;
        self.conductivity_at_temperature(to)
    }
}

#[test]
fn second_stage_rechecks_original_conductivity_support_and_ordered_time() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [SolidMaterial::new("bounded original law", 1.0)];
    let law = BoundedOriginalLaw(ConductivityTensor::new(IDENTITY).unwrap());
    let first = cooling_setup();
    let mut unsupported = first.clone();
    unsupported.faces[0] = ThermalFace::Temperature(500.0);
    let previous = [400.0];
    let solve = |second: &ThermalSetup, start: f64, dt: f64| {
        advance_nonlinear_tensor_conduction_sdirk2(
            &domain,
            &material,
            [&first, second],
            &[&law],
            &[2.0],
            &previous,
            start,
            dt,
            &controls(),
            &CancelGate::new(),
        )
    };
    assert!(solve(&first, 0.0, 1.0).is_ok());
    assert!(matches!(
        solve(&unsupported, 0.0, 1.0),
        Err(ChtError::InvalidInput { .. })
    ));
    for (start, dt) in [(1e16, 2.0), (-1.0, 1.0), (0.0, 0.0), (0.0, f64::INFINITY)] {
        assert!(matches!(
            solve(&first, start, dt),
            Err(ChtError::InvalidInput { .. })
        ));
    }
    assert_eq!(previous, [400.0]);
    assert_eq!(first, cooling_setup());
}
