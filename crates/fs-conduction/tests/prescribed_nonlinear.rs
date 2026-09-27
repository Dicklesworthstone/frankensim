//! Prescribed endpoint histories through the production nonlinear thermal step.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::material::ConductivityTable;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 83, kernel_id: 719, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        .unwrap()
}

fn boundary(mesh: &ConductionMesh, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("driven-face", |face| (face.centroid.iter().sum::<f64>() - 1.0).abs() < 1e-12,
            ThermalBc::dirichlet(temperature).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-12, ..LinearConfig::default() },
        energy_tolerance_j: 1e-8 }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e}");
}

fn nonlinear_material() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(280.0, 1.6), (360.0, 3.2)]).unwrap())
}

#[test]
fn nonlinear_driven_face_uses_endpoint_conductivity_and_full_old_storage() {
    let mesh = mesh();
    let boundary = boundary(&mesh, 330.0);
    let material = nonlinear_material();
    let source = ScalarField::Uniform(0.0);
    let old = [300.0; 4];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None };
        assert!(engine.advance_nonlinear(cx, p, None, &old, 0.25, config(), NonlinearStepConfig::default()).is_err());
        let result = engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(),
            NonlinearStepConfig::default()).unwrap();
        // Independent scalar equation: C_00=15/24=0.625, K_00=k(Tmean)/2.
        // At [310,330,330,330], Tmean=325 and k=2.5. Thus the free
        // residual is 0.625*10 + 0.25*(2.5/2)*(310-330) = 0.
        close(result.step.temperature[0], 310.0, 1e-8);
        assert_eq!(&result.step.temperature[1..], &[330.0; 3]);
        close(result.step.stored_energy_change_j, 62.5, 1e-8);
        close(result.step.dirichlet_in_w, 250.0, 1e-7);
        // The initial residual is evaluated with NEW fixed values, not at
        // the old equilibrium (which would have reported zero).
        close(result.initial_residual_j, 0.25 * (2.45 / 2.0) * 30.0, 1e-10);
        assert!(result.nonlinear_iterations > 1);
        assert!(result.residual_j <= result.threshold_j);
        assert!(result.step.energy_residual_j.abs() <= config().energy_tolerance_j);
        assert!(result.step.krylov_iterations <= config().linear.max_iterations);
        let frozen = ConductivityModel::isotropic_declared(2.0).unwrap();
        let wrong = engine.advance_prescribed(cx, ConductionProblem { material: &frozen, ..p },
            None, &old, 0.25, config()).unwrap();
        assert!((wrong.temperature[0] - result.step.temperature[0]).abs() > 1.0);
        assert_eq!(old, [300.0; 4]);
    });
}

#[test]
fn prescribed_nonlinear_path_agrees_with_linear_and_strict_limits() {
    let mesh = mesh();
    let driven = boundary(&mesh, 330.0);
    let fixed = boundary(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let nonlinear_material = nonlinear_material();
    let source = ScalarField::Uniform(0.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &driven, material: &material,
            source: &source, element_materials: None };
        let a = engine.advance_prescribed(cx, p, None, &[300.0; 4], 0.25, config()).unwrap();
        let b = engine.advance_nonlinear_prescribed(cx, p, None, &[300.0; 4], 0.25, config(),
            NonlinearStepConfig::default()).unwrap();
        for (&a, &b) in a.temperature.iter().zip(&b.step.temperature) { close(a, b, 1e-9); }
        close(a.dirichlet_in_w, b.step.dirichlet_in_w, 1e-8);
        let p = ConductionProblem { boundary: &fixed, material: &nonlinear_material, ..p };
        let old = [310.0, 300.0, 300.0, 300.0];
        let a = engine.advance_nonlinear(cx, p, None, &old, 0.25, config(), NonlinearStepConfig::default()).unwrap();
        let b = engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(), NonlinearStepConfig::default()).unwrap();
        assert_eq!(a.step.temperature, b.step.temperature);
        assert_eq!(a.initial_residual_j, b.initial_residual_j);
        assert_eq!(a.step.stored_energy_change_j, b.step.stored_energy_change_j);
    });
}

#[test]
fn nonlinear_boundary_refusals_leave_history_retryable() {
    let mesh = mesh();
    let valid = boundary(&mesh, 330.0);
    let outside = boundary(&mesh, 500.0);
    let material = nonlinear_material();
    let source = ScalarField::Uniform(0.0);
    let old = [300.0; 4];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &valid, material: &material,
            source: &source, element_materials: None };
        let one = NonlinearStepConfig { max_iterations: 1, ..NonlinearStepConfig::default() };
        assert!(matches!(engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(), one),
            Err(ConductionError::NotConverged { .. })));
        assert!(matches!(engine.advance_nonlinear_prescribed(cx, ConductionProblem { boundary: &outside, ..p },
            None, &old, 0.25, config(), NonlinearStepConfig::default()),
            Err(ConductionError::OutsideTemperatureSpan { .. })));
        let retry = engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(),
            NonlinearStepConfig::default()).unwrap();
        let fresh = engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(),
            NonlinearStepConfig::default()).unwrap();
        assert_eq!(retry.step.temperature, fresh.step.temperature);
        gate.request();
        assert!(matches!(engine.advance_nonlinear_prescribed(cx, p, None, &old, 0.25, config(),
            NonlinearStepConfig::default()), Err(ConductionError::Cancelled { .. })));
    });
    assert_eq!(old, [300.0; 4]);
}

#[test]
fn loose_nonlinear_tolerance_cannot_hide_driven_boundary_energy_defect() {
    let mesh = mesh();
    let boundary = boundary(&mesh, 330.0);
    let material = nonlinear_material();
    let source = ScalarField::Uniform(0.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None };
        let loose = NonlinearStepConfig { residual_atol_j: 1e6, ..NonlinearStepConfig::default() };
        let error = engine.advance_nonlinear_prescribed(cx, p, None, &[300.0; 4], 0.25,
            config(), loose).unwrap_err();
        assert!(error.to_string().contains("energy residual"));
    });
}
