use super::*;
use crate::adjoint::robin::RobinLinearization;
use crate::{ConductionMesh, ConductivityModel, ConductivityTable, InitialGuess,
    ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use crate::fixtures::{box_grid, on_box_face};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 59, kernel_id: 847, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    let (cells, positions) = box_grid([3,2,2], [0.1,0.04,0.03]);
    ConductionMesh::new(cells, positions).unwrap()
}
fn physical(cx: &Cx<'_>, fixed: f64, nonlinear: bool, direct: bool) -> (f64, f64) {
    let mesh = mesh();
    let high = if nonlinear {21.0} else {1.0};
    let model = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0,1.0),(450.0,high)]).unwrap());
    let source = ScalarField::Uniform(2000.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed", |f| on_box_face(f.centroid[0],0.0), ThermalBc::dirichlet(fixed).unwrap()).unwrap()
        .remainder("cooled", ThermalBc::robin(80.0,290.0).unwrap()).unwrap().finish().unwrap();
    let problem = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&model,
        element_materials:None,source:&source};
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(320.0);
    config.linear.tolerance = 1e-11;
    config.stop.residual_rtol = 1e-12;
    config.stop.step_atol = 0.0;
    let response = RobinLinearization::new(cx,problem,config,&["cooled"]).unwrap();
    let field = &response.primal().temperature;
    let saved = field.clone();
    let selected = if direct {boundary.dirichlet()[0].0} else {mesh.vertex_count()-1};
    let mut weights = vec![0.0;mesh.vertex_count()]; weights[selected] = 1.0;
    let dual = response.pullback(cx,&weights,&[0.0],&[0.0]).unwrap();
    let bars = RobinResponse::prescribed_temperature_pullback(cx,problem,None,field,
        &weights,&dual.nodal_load,&[],&[],&[],&[],1_000_000).unwrap();
    let derivative = boundary.dirichlet().iter().map(|&(v,_)| bars[v]).sum();
    for v in 0..bars.len() {
        if !boundary.dirichlet().iter().any(|&(i,_)| i == v) { assert_eq!(bars[v],0.0); }
    }
    assert_eq!(field,&saved,"no substitute primal or field mutation");
    (field[selected],derivative)
}

#[test]
fn prescribed_controls_include_full_nonlinear_lift_and_direct_fixed_goal() {
    context(&CancelGate::new_clock_free(),|cx| {
        for nonlinear in [false,true] { for direct in [false,true] {
            let (value,actual) = physical(cx,350.0,nonlinear,direct);
            let step = 0.001;
            let expected = (physical(cx,350.0+step,nonlinear,direct).0
                - physical(cx,350.0-step,nonlinear,direct).0)/(2.0*step);
            assert!((actual-expected).abs() < 2e-6*expected.abs().max(0.01),
                "nonlinear={nonlinear} direct={direct}: {actual} vs {expected}");
            if direct { assert_eq!(value,350.0); assert_eq!(actual,1.0); }
            else { assert!(actual > 0.0 && actual < 1.0); }
        } }
    });
}

#[test]
fn prescribed_controls_reject_changed_fixed_values_nonzero_fixed_duals_and_budgets() {
    let gate = CancelGate::new_clock_free();
    context(&gate,|cx| {
        let mesh = mesh();
        let model = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0,1.0),(450.0,21.0)]).unwrap());
        let source = ScalarField::Uniform(0.0);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("fixed",|f| on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(320.0).unwrap()).unwrap()
            .remainder("cooled",ThermalBc::robin(80.0,300.0).unwrap()).unwrap().finish().unwrap();
        let problem = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&model,
            element_materials:None,source:&source};
        let t = vec![320.0;mesh.vertex_count()]; let z = vec![0.0;t.len()];
        let call = |t:&[f64],lambda:&[f64],budget| RobinResponse::prescribed_temperature_pullback(
            cx,problem,None,t,&z,lambda,&[],&[],&[],&[],budget);
        assert_eq!(call(&t,&z,1_000_000).unwrap(),z);
        assert!(call(&t,&z,1).is_err());
        assert!(call(&t[..t.len()-1],&z,1_000_000).is_err());
        let v = boundary.dirichlet()[0].0;
        let mut wrong = t.clone(); wrong[v] += 1.0;
        assert!(call(&wrong,&z,1_000_000).is_err());
        let mut lambda = z.clone(); lambda[v] = 1.0;
        assert!(call(&t,&lambda,1_000_000).is_err());
        gate.request();
        assert!(matches!(call(&t,&z,1_000_000),Err(ConductionError::Cancelled {..})));
    });
}
