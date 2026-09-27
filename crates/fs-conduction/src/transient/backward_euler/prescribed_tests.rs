use super::*;
use crate::{ConductivityModel, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
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

#[test]
fn driven_face_matches_an_independent_one_dof_equation_and_reaction() {
    let mesh = mesh();
    let boundary = boundary(&mesh, 330.0);
    assert_eq!(boundary.dirichlet().len(), 3);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let old = vec![300.0; 4];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None };
        // Strict callers still refuse a changed boundary, rather than silently
        // opting into a different history model.
        assert!(engine.advance(cx, problem, None, &old, 0.25, config()).is_err());
        let result = engine.advance_prescribed(cx, problem, None, &old, 0.25, config()).unwrap();
        // On this unit right tet, V=1/6, C_00=rho_cp/24=1/2,
        // K_00=k/2=1 and K_0j=-k/6. Therefore
        // delta_0 = dt * K_00 * 30 / (C_00 + dt * K_00) = 10.
        close(result.temperature[0], 310.0, 1e-11);
        assert_eq!(&result.temperature[1..], &[330.0; 3]);
        // All FOUR nodal capacity changes count, including prescribed nodes.
        close(result.stored_energy_change_j, 0.5 * (10.0 + 3.0 * 30.0), 1e-10);
        close(result.dirichlet_in_w, 200.0, 1e-9);
        close(result.energy_residual_j, 0.0, config().energy_tolerance_j);
        assert!(result.relative_residual < config().linear.tolerance);
        let retry = engine.advance_prescribed(cx, problem, None, &old, 0.25, config()).unwrap();
        assert_eq!(result.temperature, retry.temperature);
        assert_eq!(old, vec![300.0; 4]);
    });
}

#[test]
fn successive_heating_and_cooling_endpoints_use_the_accepted_history() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let mut old = vec![300.0; 4];
        let mut expected_free = 300.0;
        let mut accumulated_energy = 0.0;
        for (end, dt) in [(330.0, 0.25), (350.0, 0.5), (290.0, 1.0)] {
            let boundary = boundary(&mesh, end);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
                source: &source, element_materials: None };
            let result = engine.advance_prescribed(cx, problem, None, &old, dt, config()).unwrap();
            expected_free = (0.5 * expected_free + dt * end) / (0.5 + dt);
            close(result.temperature[0], expected_free, 1e-10);
            assert_eq!(&result.temperature[1..], &[end; 3]);
            close(result.stored_energy_change_j, dt * result.dirichlet_in_w, 1e-8);
            accumulated_energy += result.stored_energy_change_j;
            old = result.temperature;
        }
        close(accumulated_energy, 0.5 * old.iter().map(|t| t - 300.0).sum::<f64>(), 1e-9);
    });
}

#[test]
fn all_prescribed_body_retains_the_no_free_dofs_refusal() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("all-fixed", ThermalBc::dirichlet(340.0).unwrap()).unwrap().finish().unwrap();
        let result = engine.advance_prescribed(cx, ConductionProblem { mesh: &mesh,
            boundary: &boundary, material: &material, source: &source, element_materials: None },
            None, &[300.0; 4], 2.0, config());
        assert!(matches!(result, Err(ConductionError::NoFreeDofs)));
    });
}

#[test]
fn unchanged_prescribed_values_match_the_strict_path() {
    let mesh = mesh();
    let boundary = boundary(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(5.0);
    let old = [350.0, 300.0, 300.0, 300.0];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None };
        let strict = engine.advance(cx, p, None, &old, 0.25, config()).unwrap();
        let prescribed = engine.advance_prescribed(cx, p, None, &old, 0.25, config()).unwrap();
        assert_eq!(strict.temperature, prescribed.temperature);
        assert_eq!(strict.dirichlet_in_w, prescribed.dirichlet_in_w);
        assert_eq!(strict.stored_energy_change_j, prescribed.stored_energy_change_j);
        assert_eq!(strict.krylov_iterations, prescribed.krylov_iterations);
    });
}

#[test]
fn invalid_or_cancelled_endpoint_changes_publish_no_state() {
    let mesh = mesh();
    let boundary = boundary(&mesh, 330.0);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let old = vec![300.0; 4];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None };
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(engine.advance_prescribed(cx, p, None, &old, dt, config()).is_err());
        }
        assert!(engine.advance_prescribed(cx, p, None, &old[..3], 0.25, config()).is_err());
        let mut bad = old.clone(); bad[0] = f64::NAN;
        assert!(engine.advance_prescribed(cx, p, None, &bad, 0.25, config()).is_err());
        let mut exhausted = config(); exhausted.linear.max_iterations = 0;
        assert!(engine.advance_prescribed(cx, p, None, &old, 0.25, exhausted).is_err());
        gate.request();
        assert!(matches!(engine.advance_prescribed(cx, p, None, &old, 0.25, config()),
            Err(ConductionError::Cancelled { .. })));
    });
    assert_eq!(old, vec![300.0; 4]);
}
