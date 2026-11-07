//! Physical nonlinear tensor storage through the native conduction owner.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CellSink, ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance,
    RotatedSolidConductivity, SolidMaterial, TensorConductivityLaw, TensorNonlinearConfig,
    ThermalFace, ThermalSetup, Voxel, VoxelDomain, advance_nonlinear_tensor_conduction,
    advance_tensor_conduction,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];

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

fn linear_in_temperature() -> SolidMaterial {
    SolidMaterial::new("original k(T)=T/300", 1.0)
        .with_conductivity_table(&[(200.0, 2.0 / 3.0), (600.0, 2.0)])
}

#[test]
fn nonlinear_cooling_resolves_original_integral_and_original_storage() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let setup = ThermalSetup::new(faces);
    let previous = [400.0];
    let capacity = [2.0];
    let dt = 1.0;
    let step = advance_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[&law],
        &capacity,
        &previous,
        4.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // Unit area and a half-cell distance .5:
    // Q_out = 2*integral_300^T (s/300) ds = (T²-300²)/300.
    // C*(T-400)/dt + Q_out = 0 gives this independent quadratic root.
    let b = 300.0 * capacity[0] / dt;
    let c = b * previous[0] + 300.0_f64.powi(2);
    let temperature = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
    let heat = (temperature * temperature - 300.0_f64.powi(2)) / 300.0;
    let solved = &step.solution.solution;
    close(solved.temperature[0], temperature);
    close(solved.face_temperature_k[0][0], 300.0);
    close(solved.outward_heat_w[0][0], heat);
    for face in 1..6 {
        close(solved.outward_heat_w[0][face], 0.0);
    }
    close(
        step.record.stored_energy_change_j,
        capacity[0] * (temperature - previous[0]),
    );
    close(step.record.boundary_outflow_j, dt * heat);
    close(step.record.source_j, 0.0);
    close(step.record.closure_j, 0.0);
    assert_eq!(step.record.time_s, 5.0);
    assert_eq!(step.record.iterations, solved.report.iterations);
    assert!(step.solution.newton_iterations > 0);
    assert_eq!(previous, [400.0]);
}

#[test]
fn nonlinear_finite_contact_pair_equilibrates_with_unequal_capacities() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let capacity = [2.0, 3.0];
    let previous = [400.0, 300.0];
    let dt = 2.0;
    let resistance = 0.5;
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: resistance,
    });
    let step = advance_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[&law, &law],
        &capacity,
        &previous,
        0.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // Independent scalar interface solve, using exact bulk integrals and the
    // original finite contact jump. No fitted conductivity or frozen mean.
    let temperatures = |heat: f64| {
        [
            previous[0] - dt * heat / capacity[0],
            previous[1] + dt * heat / capacity[1],
        ]
    };
    let surfaces = |heat: f64| {
        let t = temperatures(heat);
        [
            (t[0] * t[0] - 300.0 * heat).sqrt(),
            (t[1] * t[1] + 300.0 * heat).sqrt(),
        ]
    };
    let mut lo = 0.0;
    let mut hi = (previous[0] - previous[1]) / (dt * (1.0 / capacity[0] + 1.0 / capacity[1]));
    for _ in 0..80 {
        let heat = 0.5 * lo + 0.5 * hi;
        let surface = surfaces(heat);
        if surface[0] - surface[1] > resistance * heat {
            lo = heat;
        } else {
            hi = heat;
        }
    }
    let heat = 0.5 * lo + 0.5 * hi;
    let expected = temperatures(heat);
    let face = surfaces(heat);
    let solved = &step.solution.solution;
    for cell in 0..2 {
        close(solved.temperature[cell], expected[cell]);
    }
    close(solved.face_temperature_k[0][1], face[0]);
    close(solved.face_temperature_k[1][0], face[1]);
    close(solved.outward_heat_w[0][1], heat);
    close(solved.outward_heat_w[1][0], -heat);
    close(face[0] - face[1], resistance * heat);
    close(step.record.stored_energy_change_j, 0.0);
    close(step.record.source_j, 0.0);
    close(step.record.boundary_outflow_j, 0.0);
    close(step.record.closure_j, 0.0);
    assert!(solved.temperature[0] < previous[0]);
    assert!(solved.temperature[1] > previous[1]);
    assert_eq!(previous, [400.0, 300.0]);
}

#[test]
fn nonlinear_cooling_retains_physical_first_order_time_convergence() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let setup = ThermalSetup::new(faces);
    let config = controls();
    let gate = CancelGate::new();
    let march = |count: usize| {
        let dt = 1.0 / count as f64;
        let mut previous = vec![400.0];
        for index in 0..count {
            let step = advance_nonlinear_tensor_conduction(
                &domain,
                &material,
                &setup,
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
            previous = step.solution.solution.temperature;
        }
        previous[0]
    };
    // C*dT/dt = -(T²-Tb²)/300:
    // (T-Tb)/(T+Tb) decays as exp(-2*Tb*t/(300*C)).
    let ratio = (100.0 / 700.0) * (-1.0_f64).exp();
    let exact = 300.0 * (1.0 + ratio) / (1.0 - ratio);
    let coarse = (march(20) - exact).abs();
    let fine = (march(40) - exact).abs();
    assert!((1.9..2.1).contains(&(coarse / fine)), "{coarse} {fine}");
}

#[test]
fn constant_rotated_laws_reproduce_existing_transient_heat_and_storage() {
    let domain = VoxelDomain::graded_from_fn([vec![0.1, 0.3, 0.6], vec![0.2], vec![0.3]], |_| {
        Voxel::Solid(0)
    })
    .unwrap();
    let material = [SolidMaterial::new("constant owner", 1.0)];
    let tensor =
        ConductivityTensor::new([[5.0, 1.0, 0.3], [1.0, 3.0, 0.2], [0.3, 0.2, 2.0]]).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(400.0);
    faces[1] = ThermalFace::Convective {
        h: 8.0,
        ambient: 330.0,
    };
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![1.0, 0.5, 2.0];
    setup.cell_sinks.push(CellSink {
        cell: 1,
        conductance_w_k: 0.4,
        temperature_k: 310.0,
    });
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.05,
    });
    let capacity = [0.5, 2.0, 3.0];
    let previous = [350.0, 330.0, 320.0];
    let config = controls();
    let frozen = advance_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[tensor; 3],
        &capacity,
        &previous,
        2.0,
        0.7,
        &config.linear,
        &CancelGate::new(),
    )
    .unwrap();
    let laws: [&dyn TensorConductivityLaw; 3] = [&tensor; 3];
    let original = advance_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &laws,
        &capacity,
        &previous,
        2.0,
        0.7,
        &config,
        &CancelGate::new(),
    )
    .unwrap();
    for cell in 0..3 {
        close(
            original.solution.solution.temperature[cell],
            frozen.solution.temperature[cell],
        );
        for face in 0..6 {
            close(
                original.solution.solution.face_temperature_k[cell][face],
                frozen.solution.face_temperature_k[cell][face],
            );
            close(
                original.solution.solution.outward_heat_w[cell][face],
                frozen.solution.outward_heat_w[cell][face],
            );
        }
    }
    close(
        original.record.stored_energy_change_j,
        frozen.record.stored_energy_change_j,
    );
    close(original.record.source_j, frozen.record.source_j);
    close(
        original.record.boundary_outflow_j,
        frozen.record.boundary_outflow_j,
    );
    close(original.record.closure_j, 0.0);
    assert_eq!(original.record.time_s, frozen.record.time_s);
}

#[test]
fn rotated_nonlinear_insulated_heating_stores_explicit_cell_capacities_once() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 0.25, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature().with_orthotropic([6.0, 2.0, 3.0])];
    let law = RotatedSolidConductivity::new(&material[0], ROTATION).unwrap();
    let capacity = [1.0, 2.0, 4.0]; // J/K, already including actual cell volume.
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![5.0, 10.0, 20.0];
    let step = advance_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[&law, &law, &law],
        &capacity,
        &[300.0; 3],
        0.0,
        2.0,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    for cell in 0..3 {
        close(step.solution.solution.temperature[cell], 310.0);
        for &heat in &step.solution.solution.outward_heat_w[cell] {
            close(heat, 0.0);
        }
    }
    close(step.record.stored_energy_change_j, 70.0);
    close(step.record.source_j, 70.0);
    close(step.record.boundary_outflow_j, 0.0);
    close(step.record.closure_j, 0.0);
}

#[test]
fn nonlinear_prescribed_cell_reservoir_supplies_its_storage_jump_once() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature().with_orthotropic([6.0, 2.0, 3.0])];
    let law = RotatedSolidConductivity::new(&material[0], ROTATION).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature.push((0, 350.0));
    // Native prescribed-cell powers are excluded from the free-cell source ledger.
    setup.power_w = vec![1000.0];
    let step = advance_nonlinear_tensor_conduction(
        &domain,
        &material,
        &setup,
        &[&law],
        &[7.0],
        &[300.0],
        8.0,
        2.0,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close(step.solution.solution.temperature[0], 350.0);
    close(step.record.stored_energy_change_j, 350.0);
    close(step.record.source_j, 0.0);
    close(step.record.boundary_outflow_j, -350.0);
    close(step.record.closure_j, 0.0);
    assert_eq!(step.record.time_s, 10.0);
}

struct GappedSource(ConductivityTensor);
impl TensorConductivityLaw for GappedSource {
    fn conductivity_at_temperature(&self, t: f64) -> Result<ConductivityTensor, ChtError> {
        if t.is_finite() && ((250.0..=340.0).contains(&t) || (360.0..=500.0).contains(&t)) {
            Ok(self.0)
        } else {
            Err(ChtError::InvalidInput {
                field: "test.original_tensor_source",
                reason: "unsupported temperature".into(),
            })
        }
    }

    fn mean_conductivity(&self, a: f64, b: f64) -> Result<ConductivityTensor, ChtError> {
        self.conductivity_at_temperature(a)?;
        self.conductivity_at_temperature(b)?;
        if a.min(b) < 360.0 && a.max(b) > 340.0 {
            return Err(ChtError::InvalidInput {
                field: "test.original_tensor_source",
                reason: "unsupported interval".into(),
            });
        }
        Ok(self.0)
    }
}

#[test]
fn storage_time_source_support_and_solver_refusals_preserve_inputs() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = [linear_in_temperature()];
    let law = RotatedSolidConductivity::new(&material[0], IDENTITY).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let setup = ThermalSetup::new(faces);
    let previous = [400.0];
    let config = controls();
    let solve = |laws: &[&dyn TensorConductivityLaw],
                 capacity: &[f64],
                 previous: &[f64],
                 start: f64,
                 dt: f64,
                 config: &TensorNonlinearConfig,
                 gate: &CancelGate| {
        advance_nonlinear_tensor_conduction(
            &domain, &material, &setup, laws, capacity, previous, start, dt, config, gate,
        )
    };
    for capacity in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            solve(
                &[&law],
                &[capacity],
                &previous,
                0.0,
                1.0,
                &config,
                &CancelGate::new()
            )
            .is_err()
        );
    }
    for (start, dt) in [(0.0, 0.0), (-1.0, 1.0), (1e16, 1.0), (0.0, f64::INFINITY)] {
        assert!(
            solve(
                &[&law],
                &[2.0],
                &previous,
                start,
                dt,
                &config,
                &CancelGate::new()
            )
            .is_err()
        );
    }
    assert!(
        solve(
            &[],
            &[2.0],
            &previous,
            0.0,
            1.0,
            &config,
            &CancelGate::new()
        )
        .is_err()
    );
    assert!(
        solve(
            &[&law],
            &[],
            &previous,
            0.0,
            1.0,
            &config,
            &CancelGate::new()
        )
        .is_err()
    );
    assert!(
        solve(
            &[&law],
            &[2.0],
            &[0.0],
            0.0,
            1.0,
            &config,
            &CancelGate::new()
        )
        .is_err()
    );
    let no_newton = TensorNonlinearConfig {
        max_newton_iterations: 0,
        ..config
    };
    assert!(matches!(
        solve(
            &[&law],
            &[2.0],
            &previous,
            0.0,
            1.0,
            &no_newton,
            &CancelGate::new()
        ),
        Err(ChtError::SolverNotConverged { .. }),
    ));
    let gapped = GappedSource(
        ConductivityTensor::new([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]).unwrap(),
    );
    assert!(matches!(
        solve(
            &[&gapped],
            &[2.0],
            &previous,
            0.0,
            1.0,
            &config,
            &CancelGate::new()
        ),
        Err(ChtError::InvalidInput { .. }),
    ));
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        solve(&[&law], &[2.0], &previous, 0.0, 1.0, &config, &gate).unwrap_err(),
        ChtError::Cancelled,
    );
    assert_eq!(previous, [400.0]);
    assert_eq!(setup.faces[0], ThermalFace::Temperature(300.0));
}
