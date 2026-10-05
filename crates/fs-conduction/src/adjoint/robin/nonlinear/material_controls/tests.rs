use super::*;
use crate::adjoint::robin::RobinLinearization;
use crate::{ConductionMesh, ConductivityModel, ConductivityTable, ElementMaterials,
    InitialGuess, MaterialId, MaterialTable, ScalarField, SolveConfig, ThermalBc,
    ThermalBoundaryBuilder};
use crate::fixtures::{box_grid, on_box_face};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 53, kernel_id: 839, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    let (cells, positions) = box_grid([4, 2, 2], [0.1, 0.04, 0.03]);
    ConductionMesh::new(cells, positions).unwrap()
}
fn config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(320.0);
    config.linear.tolerance = 1e-11;
    config.stop.residual_rtol = 1e-12;
    config.stop.step_atol = 0.0;
    config
}
fn physical(cx: &Cx<'_>, scales: [f64; 2], nonlinear: bool, differentiated: bool) -> (f64, Vec<f64>) {
    let mesh = mesh();
    let model = |i: usize| {
        let values = if nonlinear {
            if i == 0 { [1.0,21.0] } else { [6.0,14.0] }
        } else if i == 0 { [9.0,9.0] } else { [13.0,13.0] };
        ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![
            (250.0,values[0]*scales[i]),(450.0,values[1]*scales[i])]).unwrap())
    };
    // Assignment ownership, not the unusable fallback, selects every tensor.
    let fallback = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![
        (0.0,1.0),(1.0,2.0)]).unwrap());
    let materials = ElementMaterials::new(MaterialTable::new([
        (MaterialId(7),model(0)),(MaterialId(8),model(1))]).unwrap(),
        (0..mesh.element_count()).map(|e| MaterialId(7+(e%2) as u32)).collect()).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |f| on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(350.0).unwrap()).unwrap()
        .region("cooled", |f| on_box_face(f.centroid[0],0.1),ThermalBc::robin(80.0,290.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let source = ScalarField::Uniform(2000.0);
    let problem = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&fallback,
        element_materials:Some(&materials),source:&source};
    let linearization = RobinLinearization::new(cx,problem,config(),&["cooled"]).unwrap();
    let field = &linearization.primal().temperature;
    let value = field[mesh.vertex_count()-1];
    if !differentiated { return (value,Vec::new()); }
    let mut weights = vec![0.0;mesh.vertex_count()];
    weights[mesh.vertex_count()-1]=1.0;
    let gradient = linearization.pullback(cx,&weights,&[0.0],&[0.0]).unwrap();
    let entries = RobinResponse::conductivity_scale_pullback(cx,problem,field,
        &gradient.nodal_load,mesh.element_count()).unwrap();
    let mut grouped = vec![0.0;2];
    for (e,value) in entries.into_iter().enumerate() { grouped[e%2]+=value; }
    assert_eq!(field,&linearization.primal().temperature,"no replacement primal");
    (value,grouped)
}

#[test]
fn conductivity_scales_match_heterogeneous_linear_and_nonlinear_physical_resolves() {
    context(&CancelGate::new_clock_free(),|cx| {
        for nonlinear in [false,true] {
            let scales=[0.8,1.3];
            let (_,actual)=physical(cx,scales,nonlinear,true);
            let epsilon=2e-4_f64;
            for i in 0..2 {
                let mut low=scales;low[i]*=(-epsilon).exp();
                let mut high=scales;high[i]*=epsilon.exp();
                let expected=(physical(cx,high,nonlinear,false).0-physical(cx,low,nonlinear,false).0)/(2.0*epsilon);
                assert!(expected.abs()>0.01,"fixture must observe each independent material");
                assert!((actual[i]-expected).abs()<3e-5*expected.abs().max(1.0),
                    "nonlinear={nonlinear} region={i}: {} vs {expected}",actual[i]);
            }
        }
    });
}

#[test]
fn conductivity_controls_refuse_wrong_fields_fixed_duals_budgets_and_kinks() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();
        let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![
            (250.0,1.0),(450.0,21.0)]).unwrap());
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("fixed",|f| on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(320.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
            element_materials:None,source:&source};
        let t=vec![320.0;mesh.vertex_count()];let z=vec![0.0;t.len()];let budget=mesh.element_count();
        assert_eq!(RobinResponse::conductivity_scale_pullback(cx,problem,&t,&z,budget).unwrap(),vec![0.0;budget]);
        assert!(RobinResponse::conductivity_scale_pullback(cx,problem,&t,&z,budget-1).is_err());
        assert!(RobinResponse::conductivity_scale_pullback(cx,problem,&t[..t.len()-1],&z,budget).is_err());
        let fixed=boundary.dirichlet()[0].0;
        let mut wrong=z.clone();wrong[fixed]=1.0;
        assert!(RobinResponse::conductivity_scale_pullback(cx,problem,&t,&wrong,budget).is_err());
        let mut wrong_t=t.clone();wrong_t[fixed]=321.0;
        assert!(RobinResponse::conductivity_scale_pullback(cx,problem,&wrong_t,&z,budget).is_err());
        wrong[fixed]=0.0;wrong[t.len()-1]=f64::NAN;
        assert!(RobinResponse::conductivity_scale_pullback(cx,problem,&t,&wrong,budget).is_err());
        let kink=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![
            (250.0,1.0),(320.0,4.0),(450.0,21.0)]).unwrap());
        assert!(RobinResponse::conductivity_scale_pullback(cx,ConductionProblem {material:&kink,..problem},&t,&z,budget).is_err());
    });
    let gate=CancelGate::new_clock_free();gate.request();
    context(&gate,|cx| {
        let mesh=mesh();
        let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0,1.0),(450.0,21.0)]).unwrap());
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .remainder("cooled",ThermalBc::robin(80.0,300.0).unwrap()).unwrap().finish().unwrap();
        let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,element_materials:None,source:&source};
        assert!(matches!(RobinResponse::conductivity_scale_pullback(cx,problem,&vec![300.0;mesh.vertex_count()],
            &vec![0.0;mesh.vertex_count()],mesh.element_count()),Err(ConductionError::Cancelled {..})));
    });
}
