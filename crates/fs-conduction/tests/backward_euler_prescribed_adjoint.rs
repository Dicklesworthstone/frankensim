//! Prescribed-boundary launch states must retain storage in their adjoints.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, StepConfig, StepLinearization,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 71, kernel_id: 824, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn config() -> StepConfig {
    StepConfig {
        linear: LinearConfig { tolerance: 1e-11, max_iterations: 2000, restart: 20 },
        energy_tolerance_j: 1e-7,
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual:e} differs from {expected:e}");
}

// Controls are the initial free temperature, volumetric source multiplier,
// and common capacity multiplier. The prescribed endpoint stays at 340 K.
fn case<R>(cx: &Cx<'_>, control: [f64; 3], nonlinear: bool,
    f: impl FnOnce(&StepLinearization<'_>) -> R) -> R {
    let (complex, positions) = box_grid([2, 1, 1], [0.1, 0.04, 0.03]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("heater", |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).unwrap()).unwrap()
        .region("air", |face| on_box_face(face.centroid[0], 0.1),
            ThermalBc::robin(40.0, 290.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = if nonlinear {
        ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 1.0), (450.0, 21.0)]).unwrap())
    } else { ConductivityModel::isotropic_declared(10.0).unwrap() };
    let source = ScalarField::Uniform(2000.0 * control[1]);
    let old: Vec<_> = mesh.positions().iter().map(|p| {
        // Fixed old nodes are independent history, not a prescribed control.
        if on_box_face(p[0], 0.0) { 300.0 } else { 300.0 + control[0] }
    }).collect();
    let saved = old.clone();
    let engine = BackwardEuler::uniform(cx, &mesh,
        VolumetricHeatCapacity::declared(1000.0 * control[2]).unwrap()).unwrap();
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let policy = nonlinear.then_some(NonlinearStepConfig::default());
    let direct = match policy {
        Some(policy) => engine.advance_nonlinear_prescribed(
            cx, problem, None, &old, 2.0, config(), policy).unwrap().step,
        None => engine.advance_prescribed(cx, problem, None, &old, 2.0, config()).unwrap(),
    };
    let step = engine.linearize_prescribed_step(cx, problem, None, &old, 2.0,
        config(), policy, &["air"]).unwrap();
    assert_eq!(old, saved);
    assert_eq!(step.primal().temperature, direct.temperature);
    assert_eq!(step.primal().stored_energy_change_j.to_bits(), direct.stored_energy_change_j.to_bits());
    assert_eq!(step.primal().dirichlet_in_w.to_bits(), direct.dirichlet_in_w.to_bits());
    assert!(step.primal().energy_residual_j.abs() <= config().energy_tolerance_j);
    f(&step)
}
fn mean(step: &StepLinearization<'_>) -> f64 {
    step.temperature().iter().sum::<f64>() / step.temperature().len() as f64
}

#[test]
fn prescribed_launch_history_source_and_capacity_match_finite_differences() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        for nonlinear in [false, true] {
            let nominal = [0.0, 1.0, 1.0];
            let gradient = case(cx, nominal, nonlinear, |step| {
                let weights = vec![1.0 / step.temperature().len() as f64; step.temperature().len()];
                let g = step.pullback(cx, &weights, &[0.0], &[0.0]).unwrap();
                assert!(g.relative_residual < config().linear.tolerance);
                [step.previous_temperature_pullback(cx, &g.nodal_load).unwrap().iter().sum(),
                    step.source_multiplier_pullback(cx, &g.nodal_load).unwrap(),
                    step.capacity_multiplier_pullback(cx, &g.nodal_load).unwrap()]
            });
            for axis in 0..3 {
                let h = 1e-4;
                let mut plus = nominal; plus[axis] += h;
                let mut minus = nominal; minus[axis] -= h;
                let expected = (case(cx, plus, nonlinear, mean)
                    - case(cx, minus, nonlinear, mean)) / (2.0 * h);
                close(gradient[axis], expected, 5e-5);
            }
            assert!(gradient[0] > 0.0 && gradient[0] < 1.0);
            assert!(gradient[1] > 0.0);
            assert!(gradient[2].abs() > 1e-3, "capacity must remain in the endpoint Jacobian");
        }
    });
}

#[test]
fn prescribed_launch_does_not_change_legacy_admission_or_ignore_cancellation() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let (complex, positions) = box_grid([1, 1, 1], [0.1, 0.04, 0.03]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("heater", |face| on_box_face(face.centroid[0], 0.0),
                ThermalBc::dirichlet(340.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let old = vec![300.0; mesh.vertex_count()];
        let engine = BackwardEuler::uniform(cx, &mesh,
            VolumetricHeatCapacity::declared(1000.0).unwrap()).unwrap();
        assert!(engine.linearize_step(cx, problem, None, &old, 2.0, config(), None, &[]).is_err());
        let gate = CancelGate::new_clock_free();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(engine.linearize_prescribed_step(
                cancelled, problem, None, &old, 2.0, config(), None, &[]).is_err());
        });
        for dt in [0.0, -1.0, f64::NAN] {
            assert!(engine.linearize_prescribed_step(cx, problem, None, &old, dt,
                config(), None, &[]).is_err());
        }
        assert_eq!(old, vec![300.0; mesh.vertex_count()]);
    });
}
