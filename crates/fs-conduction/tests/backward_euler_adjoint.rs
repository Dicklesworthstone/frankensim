use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig, StepLinearization};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 61, kernel_id: 820, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R { with_gate(&CancelGate::new_clock_free(), f) }
fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, max_iterations: 2000, restart: 20 },
        energy_tolerance_j: 1e-8 }
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol * b.abs().max(1.0), "{a:e} != {b:e}"); }
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }

// controls: uniform initial-temperature shift, source multiplier, capacity
// multiplier, hot reference shift, hot log(h) change.
fn case<R>(cx: &Cx<'_>, control: [f64; 5], nonlinear: bool,
    f: impl FnOnce(&StepLinearization<'_>) -> R) -> R {
    let (complex, positions) = box_grid([2, 1, 1], [0.1, 0.04, 0.03]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::robin(40.0 * control[4].exp(), 320.0 + control[3]).unwrap()).unwrap()
        .region("cold", |face| on_box_face(face.centroid[0], 0.1), ThermalBc::robin(80.0, 290.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = if nonlinear {
        ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0, 1.0), (450.0, 21.0)]).unwrap())
    } else { ConductivityModel::isotropic_declared(10.0).unwrap() };
    let source = ScalarField::Uniform(2000.0 * control[1]);
    let old: Vec<_> = mesh.positions().iter().map(|p| 300.0 + 100.0 * p[0] + control[0]).collect();
    let engine = BackwardEuler::uniform(cx, &mesh,
        VolumetricHeatCapacity::declared(1000.0 * control[2]).unwrap()).unwrap();
    let step = engine.linearize_step(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source }, None, &old, 2.0,
        config(), nonlinear.then_some(NonlinearStepConfig::default()), &["hot", "cold"]).unwrap();
    f(&step)
}
fn mean(step: &StepLinearization<'_>) -> f64 {
    step.primal().temperature.iter().sum::<f64>() / step.primal().temperature.len() as f64
}

#[test]
fn history_source_and_capacity_pullbacks_match_perturbed_nonlinear_endpoints() {
    with_cx(|cx| for nonlinear in [false, true] {
        let nominal = [0.0, 1.0, 1.0, 0.0, 0.0];
        let derivatives = case(cx, nominal, nonlinear, |step| {
            let weights = vec![1.0 / step.temperature().len() as f64; step.temperature().len()];
            let g = step.pullback(cx, &weights, &[0.0;2], &[0.0;2]).unwrap();
            [step.previous_temperature_pullback(cx, &g.nodal_load).unwrap().iter().sum(),
                step.source_multiplier_pullback(cx, &g.nodal_load).unwrap(),
                step.capacity_multiplier_pullback(cx, &g.nodal_load).unwrap(),
                g.references[0], g.log_htc[0]]
        });
        for axis in 0..5 {
            let eps = 1e-4;
            let mut plus = nominal; plus[axis] += eps;
            let mut minus = nominal; minus[axis] -= eps;
            let fd = (case(cx, plus, nonlinear, mean) - case(cx, minus, nonlinear, mean)) / (2.0 * eps);
            close(derivatives[axis], fd, 3e-5);
        }
        assert!(derivatives[0] > 0.0 && derivatives[0] < 1.0,
            "transient history must neither disappear nor be an identity");
    });
}

#[test]
fn endpoint_transpose_identity_includes_direct_wall_heat_terms() {
    with_cx(|cx| case(cx, [0.0,1.0,1.0,0.0,0.0], true, |step| {
        let mut d = step.zero_direction();
        d.references_k = vec![0.3,-0.2]; d.log_htc = vec![-0.1,0.4];
        for (i,v) in d.nodal_load_w.iter_mut().enumerate() { *v = (i%3) as f64 * 0.001; }
        let weights: Vec<_> = (0..d.nodal_load_w.len()).map(|i| (i%5) as f64 * 0.03).collect();
        let tangent = step.apply(cx,&d).unwrap();
        let gradient = step.pullback(cx,&weights,&[0.4,-0.2],&[0.03,-0.01]).unwrap();
        let lhs = dot(&weights,&tangent.temperature_k) + dot(&[0.4,-0.2],&tangent.mean_wall_temperatures_k)
            + dot(&[0.03,-0.01],&tangent.heat_rates_w);
        close(lhs, dot(&gradient.references,&d.references_k) + dot(&gradient.log_htc,&d.log_htc)
            + dot(&gradient.nodal_load,&d.nodal_load_w), 1e-8);
    }));
}

#[test]
fn uniform_insulated_heating_has_exact_source_capacity_and_history_sensitivities() {
    with_cx(|cx| {
        let (complex, positions) = box_grid([1,1,1],[1.0,1.0,1.0]);
        let mesh = ConductionMesh::new(complex,positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(10.0);
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let step = engine.linearize_step(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,
            material:&material,element_materials:None,source:&source},None,&vec![300.0;mesh.vertex_count()],
            2.0,config(),None,&[]).unwrap();
        for &t in &step.primal().temperature { close(t,304.0,1e-12); }
        let w = vec![1.0/mesh.vertex_count() as f64;mesh.vertex_count()];
        let g = step.pullback(cx,&w,&[],&[]).unwrap();
        close(step.source_multiplier_pullback(cx,&g.nodal_load).unwrap(),4.0,1e-10);
        close(step.capacity_multiplier_pullback(cx,&g.nodal_load).unwrap(),-4.0,1e-10);
        close(step.previous_temperature_pullback(cx,&g.nodal_load).unwrap().iter().sum(),1.0,1e-10);
    });
}

#[test]
fn cancelled_and_malformed_pullbacks_publish_no_partial_vector() {
    with_cx(|cx| case(cx,[0.0,1.0,1.0,0.0,0.0],false,|step| {
        assert!(step.previous_temperature_pullback(cx,&[1.0]).is_err());
        assert!(step.capacity_multiplier_pullback(cx,&vec![f64::NAN;step.temperature().len()]).is_err());
        let gate = CancelGate::new_clock_free(); gate.request();
        with_gate(&gate, |blocked| {
            let zero = vec![0.0;step.temperature().len()];
            assert!(step.previous_temperature_pullback(blocked,&zero).is_err());
            assert!(step.pullback(blocked,&zero,&[0.0;2],&[0.0;2]).is_err());
        });
    }));
}
