use super::*;
use crate::{ConductivityModel, ScalarField, ThermalBc, ThermalBoundaryBuilder};
use crate::fixtures::{box_grid, on_box_face};
use crate::transient::{march, TransientConfig, TransientProblem};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 79, kernel_id: 719, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, ..LinearConfig::default() }, energy_tolerance_j: 1e-7 }
}
fn mesh() -> ConductionMesh {
    let (complex, points) = box_grid([2,2,2], [0.1,0.1,0.1]);
    ConductionMesh::new(complex, points).unwrap()
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }

#[test]
fn insulated_heating_is_exact_and_repeated_trials_do_not_advance_history() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(2000.0);
    let old = vec![300.0; mesh.vertex_count()];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let step = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(2e6).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material, source: &source, element_materials: None };
        let a = step.advance(cx, p, None, &old, 2.0, config()).unwrap();
        let b = step.advance(cx, p, None, &old, 2.0, config()).unwrap();
        assert_eq!(a.temperature, b.temperature);
        assert!(old.iter().all(|&t| t == 300.0));
        for &t in &a.temperature { close(t, 300.002, 1e-10); }
        close(a.source_w, 2.0, 1e-10);
        close(a.stored_energy_change_j, 4.0, 1e-7);
        assert!(a.energy_residual_j.abs() <= config().energy_tolerance_j);
        let next = step.advance(cx, p, None, &a.temperature, 2.0, config()).unwrap();
        for &t in &next.temperature { close(t, 300.004, 1e-10); }
    });
}

#[test]
fn one_robin_step_matches_the_existing_transient_march() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("cooling", ThermalBc::robin(80.0, 300.0).unwrap()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let old = vec![350.0; mesh.vertex_count()];
    let capacity = VolumetricHeatCapacity::declared(2e6).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, capacity).unwrap();
        let a = engine.advance(cx, ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source }, None, &old, 3.0, config()).unwrap();
        let b = march(cx, TransientProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, capacity }, &TransientConfig::backward_euler(3.0, config().linear).unwrap(), &old, 1).unwrap();
        for (a,b) in a.temperature.iter().zip(&b.temperature) { close(*a,*b,1e-7); }
        assert!(a.stored_energy_change_j < 0.0);
        close(a.stored_energy_change_j, -3.0*a.robin_out_w, 1e-7);
        close(a.robin_fluxes.iter().map(|f| f.heat_rate_w).sum(), a.robin_out_w, 1e-10);
    });
}

#[test]
fn per_element_capacities_keep_disconnected_solids_distinct() {
    let points = vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.],
        [2.,0.,0.], [3.,0.,0.], [2.,1.,0.], [2.,0.,1.]];
    let mesh = ConductionMesh::new(TetComplex::from_tets(8, vec![[0,1,2,3], [4,5,6,7]]), points).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(1.0).unwrap();
    let source = ScalarField::Uniform(2.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let capacities = [10.0,20.0].map(|v| VolumetricHeatCapacity::declared(v).unwrap());
        let engine = BackwardEuler::per_element(cx, &mesh, &capacities).unwrap();
        let result = engine.advance(cx, ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source }, None, &[300.0;8], 0.5, config()).unwrap();
        for &t in &result.temperature[..4] { close(t,300.1,1e-10); }
        for &t in &result.temperature[4..] { close(t,300.05,1e-10); }
        close(result.stored_energy_change_j, 1.0/3.0, 1e-10);
        assert!(BackwardEuler::per_element(cx, &mesh, &capacities[..1]).is_err());
    });
}

#[test]
fn dirichlet_reaction_includes_capacity_and_does_not_add_an_absolute_lift() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed", |f| on_box_face(f.centroid[0],0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(1000.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(1e4).unwrap()).unwrap();
        let p = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source };
        let result = engine.advance(cx, p, None, &vec![300.0;mesh.vertex_count()], 1.0, config()).unwrap();
        for &(v,t) in boundary.dirichlet() { assert_eq!(result.temperature[v],t); }
        close(result.stored_energy_change_j, result.source_w+result.dirichlet_in_w, 1e-7);
        let mut wrong = vec![300.0;mesh.vertex_count()]; wrong[boundary.dirichlet()[0].0]=301.0;
        assert!(engine.advance(cx,p,None,&wrong,1.0,config()).is_err());
    });
}

#[test]
fn malformed_and_cancelled_steps_return_no_field() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(1.0).unwrap();
    let source = ScalarField::Uniform(1.0);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        for dt in [0.0,-1.0,f64::NAN,f64::INFINITY] {
            assert!(engine.advance(cx,p,None,&vec![300.0;mesh.vertex_count()],dt,config()).is_err());
        }
        assert!(engine.advance(cx,p,None,&[],1.0,config()).is_err());
        let mut exhausted = config(); exhausted.linear.max_iterations=0;
        assert!(engine.advance(cx,p,None,&vec![300.0;mesh.vertex_count()],1.0,exhausted).is_err());
        gate.request();
        assert!(matches!(engine.advance(cx,p,None,&vec![300.0;mesh.vertex_count()],1.0,config()), Err(ConductionError::Cancelled{..})));
    });
}

#[test]
fn checked_correction_keeps_a_hard_iteration_cap_and_retryable_history() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("cooling", ThermalBc::robin(80.0, 300.0).unwrap()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(2000.0);
    let mut old = vec![300.0; mesh.vertex_count()]; old[0] = 350.0;
    let saved = old.clone();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(2e6).unwrap()).unwrap();
        let p = ConductionProblem {mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None};
        let mut policy = config(); policy.linear.tolerance = 1e-12;
        let mut short = policy; short.linear.max_iterations = 1;
        match engine.advance(cx, p, None, &old, 3.0, short) {
            Err(ConductionError::LinearSolveFailed {krylov_iterations, true_relative_residual, tolerance, ..}) => {
                assert_eq!(krylov_iterations, 1);
                assert_eq!(tolerance, policy.linear.tolerance);
                assert!(true_relative_residual >= tolerance);
            }
            other => panic!("expected the unchanged work cap to refuse: {other:?}"),
        }
        assert_eq!(old, saved);
        let retry = engine.advance(cx, p, None, &old, 3.0, policy).unwrap();
        let fresh = engine.advance(cx, p, None, &old, 3.0, policy).unwrap();
        assert_eq!(retry.temperature, fresh.temperature);
        assert_eq!(retry.krylov_iterations, fresh.krylov_iterations);
        assert!(retry.krylov_iterations <= policy.linear.max_iterations);
        assert!(retry.relative_residual < policy.linear.tolerance);
        assert!(retry.energy_residual_j.abs() <= policy.energy_tolerance_j);
        assert_eq!(old, saved);
    });
}
